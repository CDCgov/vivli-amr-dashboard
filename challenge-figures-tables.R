# pkgs
library(tidyverse)
library(glue)

# load in data
vivli_summ <- readr::read_csv('data/vivli_summ.csv')
vivli_select <- readr::read_csv('data/vivli_select.csv')
species_by_drug <- readr::read_csv('data/species_by_drug.csv')
cutoffs_formatted <- readr::read_csv('data/cutoffs_formatted.csv')

# tables for Dallas
species_by_drug %>%
  filter(!is.na(prop_resistant)) %>%
  mutate(table_val = glue("{round(prop_resistant * 100, 2)}% ({round(prop_resistant * n_tested)} / {n_tested})")) %>%
  select(drug, Species, table_val) %>%
  pivot_wider(names_from = drug, values_from = table_val) %>%
  write_csv("output/prop_resistant.csv")

species_by_drug %>%
  filter(!is.na(prop_nonWT)) %>%
  mutate(table_val = glue("{round(prop_nonWT * 100, 2)}% ({round(prop_nonWT * n_tested)} / {n_tested})")) %>%
  select(drug, Species, table_val) %>%
  pivot_wider(names_from = drug, values_from = table_val) %>%
  write_csv("output/prop_nonWT.csv")

species_by_drug %>%
  filter(is.na(prop_nonWT) & is.na(prop_resistant)) %>%
  mutate(table_val = modal_mic) %>%
  select(drug, Species, table_val) %>%
  pivot_wider(names_from = drug, values_from = table_val) %>%
  write_csv("output/modal_mic.csv")

# figures by region, 