library(foundry)
library(tidyverse)

# import data from dcipher ----
vivli <- datasets.read_table("vivli_sentry_data_july2023_fungi_cdc")
ecvs_clsi <- datasets.read_table("ecvs_clsi_breakpoints___table1")
whoregions <- datasets.read_table("ecvs_clsi_breakpoints___whoregions")

# clean up species ----
vivli <-
  vivli %>%
  mutate(Species = case_when(grepl("Cryptococcus neoformans", Species) ~ "Cryptococcus neoformans", 
                             TRUE ~ Species))
ecvs_clsi <-
  ecvs_clsi %>%
  mutate(Species = case_when(grepl("Cryptococcus neoformans", Species) ~ "Cryptococcus neoformans", 
                             TRUE ~ Species))

# filter to species with more than 10 isolates ----
species_for_tables <-
  vivli %>%
  count(Species) %>%
  filter(n >= 10)

# filter to select species and join with cutoffs/ecvs ----
vivli_long <-
  vivli %>%
  select(-ends_with("CLSI_CLSI")) %>%
  pivot_longer(Amphotericin_B_CLSI:Voriconazole_CLSI, 
               names_to = "drug", values_to = "mic") %>%
  mutate(drug = gsub("_CLSI", "", drug)) %>%
  mutate(drug = gsub("_", " ", drug)) %>%
  left_join(whoregions %>% select(-Year), by = c("Country" = "Entity"))

vivli_long <-
  vivli_long %>%
  left_join(ecvs_clsi, by = c("Species", "drug" = "Antifungal")) %>%
  mutate(resistance_status = case_when(
    !is.na(Resistance_lower) & mic >= Resistance_lower ~ "Resistant", 
    !is.na(SDD_upper) & mic < Resistance_lower & mic >= SDD_upper ~ "SDD", 
    !is.na(Intermediate_upper) & !is.na(SDD_upper) & mic < SDD_upper & mic >= Intermediate_upper ~ "Intermediate", 
    !is.na(Intermediate_upper) & is.na(SDD_upper) & mic < Resistance_lower & mic >= Intermediate_upper ~ "Intermediate", 
    !is.na(Sensitive_upper) & mic <= Sensitive_upper ~ "Sensitive", 
    TRUE ~ NA_character_),
    ecv_status = ifelse(mic <= `_ECV`, "WT", "non-WT")) %>%
  select(-any_of(colnames(ecvs_clsi)[colnames(ecvs_clsi) != "Species"])) 

vivli_long %>%
  summarize(n= n(), .by = c(mic, Species, drug, resistance_status, ecv_status))

vivli_select <-
  vivli_long %>%
  filter(Species %in% species_for_tables$Species)

# get modal mics ----
species_by_drug <-
  vivli_select %>%
  group_by(drug, Species) %>%
  summarize(n_tested = sum(!is.na(mic)), 
            n_isos = n(), 
            prop_nonWT = sum(ecv_status %in% "non-WT")/sum(!is.na(ecv_status)), 
            prop_resistant = sum(resistance_status %in% "Resistant") / sum(!is.na(resistance_status)))

multiple_modes <- 
  vivli_select %>%
  filter(!is.na(mic)) %>%
  group_by(drug, Species) %>%
  count(mic) %>%
  right_join(species_by_drug %>% select(drug, Species, n_tested)) %>%
  mutate(prop_tested = n/n_tested) %>%
  select(-n_tested) %>%
  slice_max(n) %>%
  ungroup() %>%  
  mutate(nmodes = n(), .by = c(Species, drug)) %>%
  filter(nmodes > 1)

write_csv(multiple_modes, "output/multiple_modal_mics.csv")

modal_mics <- 
  vivli_select %>%
  filter(!is.na(mic)) %>%
  group_by(drug, Species) %>%
  count(mic) %>%
  right_join(species_by_drug %>% select(drug, Species, n_tested)) %>%
  mutate(prop_tested = n/n_tested) %>%
  select(-n_tested) %>%
  slice_max(tibble(n, mic))


species_by_drug <- 
  species_by_drug %>%
  left_join(modal_mics  %>% 
              select(drug, Species, modal_mic = mic))


# formatted cutoff table by species x drug ----
cutoffs_formatted <-
  ecvs_clsi %>%
  select(-ends_with("upper"), -ends_with("lower"), -`_ECV`) %>%
  pivot_longer(-c(Species, Antifungal)) %>%
  group_by(Species, Antifungal) %>%
  filter(!is.na(value)) %>%
  summarize(clsi = glue::glue_collapse(glue::glue("{name}: {value}"), sep = " <br>")) %>%
  right_join(ecvs_clsi %>% select(Species, Antifungal, ecv = `_ECV`)) 

species_by_drug <-
  species_by_drug %>%
  left_join(cutoffs_formatted, by = c("Species", "drug" = "Antifungal")) 


class_checks <-
  vivli_long %>%
  summarize(n= n(), .by = c(mic, Species, drug, resistance_status, ecv_status)) %>%
  left_join(cutoffs_formatted, by = c("Species", "drug" = "Antifungal")) %>%
  filter(!is.na(clsi) | !is.na(ecv))
write_csv(class_checks, "double_check_mics.csv")

# summary statistics by different covariates ----
vivli_summ <- 
  vivli_select %>% 
  mutate(Year = as.character(Year)) %>%
  pivot_longer(c(Gender, Age_Group, WHO_region, Speciality, Source, Country, Year)) %>%
  group_by(Species, drug, name, value) %>%
  summarize(prop_resistant = sum(resistance_status %in% "Resistant")/sum(!is.na(resistance_status)),
            n_tested = sum(!is.na(mic)), 
            n_isos = n(), 
            prop_nonWT = sum(ecv_status %in% "non-WT")/sum(!is.na(ecv_status)),
            n_resistant = sum(resistance_status %in% "Resistant"), 
            tested = sum(!is.na(mic))) %>%
  mutate(
    tested = glue::glue("{n_tested}/{n_isos} ({round(n_tested/n_isos*100, 2)}%)")) %>%
  select(-n_isos)

vivli_mics <-
  vivli_select %>% 
  mutate(Year = as.character(Year)) %>%
  pivot_longer(c(Gender, Age_Group, Speciality, Source, WHO_region, Country, Year)) %>%
  group_by(Species, drug, name, value) %>%
  filter(!is.na(mic)) %>%
  count(mic) %>%
  right_join(species_by_drug %>% select(drug, Species, n_tested)) %>%
  mutate(prop_tested = n/n_tested) %>%
  select(-n_tested) %>%
  slice_max(tibble(n, mic))


vivli_mic_dists <-
  vivli_select %>% 
  mutate(Year = as.character(Year)) %>%
  pivot_longer(c(Gender, Age_Group, WHO_region, Year)) %>%
  filter(!is.na(mic)) %>%
  group_by(Species, drug, name, value, mic) %>%
  summarize(n = n(), resistance_status = resistance_status[1], 
            ecv_status = ecv_status[1]) %>%
  group_by(Species, drug, name, value) %>%
  mutate(prop = n/sum(n)) %>%
  right_join(species_by_drug %>% select(drug, Species, n_tested)) %>%
  mutate(prop_tested = n/n_tested) 

vivli_summ <- 
  vivli_summ %>%
  left_join(vivli_mics %>% select(drug, Species, name, value, modal_mic = mic))


# write out datasets for dashboard ----
readr::write_csv(vivli_summ, 'data/vivli_summ.csv')
readr::write_csv(vivli_select, 'data/vivli_select.csv')
readr::write_csv(species_by_drug, 'data/species_by_drug.csv')
readr::write_csv(cutoffs_formatted, 'data/cutoffs_formatted.csv')
readr::write_csv(vivli_mic_dists, 'data/vivli_mic_dists.csv')
