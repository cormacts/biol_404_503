#############
# BIOL 403 503 Project
# Script: 01_traits_v0
# By: Lydia, Elias, Cormac, Annie
# 2024-03-18
# This script will organize functional traits based on the data from Weigel (2022) paper
#############

library(dplyr) 

## Importing Data ####
traits <- read.csv("data/raw/msystems.01422-21-s0005.csv")

## Condensing categories of traits into single columns ####

traits <- traits %>%
  mutate_at(vars(starts_with("Dissimilatory_nitrate_reduction_")), ~coalesce(., "NA")) %>%
  mutate(dissimilatory_nitrate_reduction = if_else(rowSums(select(., starts_with("Dissimilatory_nitrate_reduction_")) == "Y", na.rm = TRUE) > 0, "Y", ""))

traits <- traits %>%
  mutate_at(vars(starts_with("Assimilatory_nitrate_reduction_")), ~coalesce(., "NA")) %>%
  mutate(assimilatory_nitrate_reduction = if_else(rowSums(select(., starts_with("Assimilatory_nitrate_reduction_")) == "Y", na.rm = TRUE) > 0, "Y", ""))

traits <- traits %>%
  mutate_at(vars(starts_with("Denitrification")), ~coalesce(., "NA")) %>%
  mutate(denitrification = if_else(rowSums(select(., starts_with("Denitrification")) == "Y", na.rm = TRUE) > 0, "Y", ""))

traits <- traits %>%
  mutate_at(vars(starts_with("Nitrogen_fixation")), ~coalesce(., "NA")) %>%
  mutate(nitrogen_fixation = if_else(rowSums(select(., starts_with("Nitrogen_fixation")) == "Y", na.rm = TRUE) > 0, "Y", ""))

traits <- traits %>%
  mutate_at(vars(starts_with("Nitrification")), ~coalesce(., "NA")) %>%
  mutate(nitrification = if_else(rowSums(select(., starts_with("Nitrification")) == "Y", na.rm = TRUE) > 0, "Y", ""))

## Condensing categories into one single nitrogen-related column ####

traits <- traits %>%
  mutate_at(vars(c("dissimilatory_nitrate_reduction","assimilatory_nitrate_reduction","denitrification","nitrogen_fixation","nitrification")), ~coalesce(., "NA")) %>%
  mutate(nitrogen_cycling = if_else(rowSums(select(., c("dissimilatory_nitrate_reduction","assimilatory_nitrate_reduction","denitrification","nitrogen_fixation","nitrification")) == "Y", na.rm = TRUE) > 0, "Y", ""))

## Creating a new condensed dataframe, with only one variable for nitrogen cycling ####

condensed_traits <- traits[1:66,c("MAG_name", "phylum","class","order","family","genus","nitrogen_cycling")]

## Finding lowest taxanomic data ####
# Creating function to remove blank strings and replace with NAs
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}

replace_na_with_n <- function(x) {
  x[is.na(x)] <- "N"
  return(x)
}

# Applying the above function to our dataset
condensed_traits <- condensed_traits %>%
  mutate(across(everything(), replace_empty_with_na))

trait_data <- condensed_traits %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(genus), genus,
                       if_else(!is.na(family), family,
                       if_else(!is.na(order), order,
                       if_else(!is.na(class), class,
                       if_else(!is.na(phylum), phylum, NA_character_))))))

trait_data <- replace_na_with_n(trait_data)




#### Attempts to propogate data ####

# columns_to_process <- c("phylum", "class", "order", "family", "genus")
# 
# tax_data <- condensed_traits %>%
#   mutate(ID = MAG_name) %>%
#   select(all_of(columns_to_process)) %>%
#   t() %>%
#   na.locf() %>%
#   t()
# 
# tax_data <- as.data.frame(tax_data)
# 
# total_data <- bind_rows(tax_data, trait_data) %>%
#   mutate(across(everything(), ~coalesce(., first(.)))) %>%
#   summarize(across(everything(), ~.))

# columns_to_process <- c("phylum", "class", "order", "family", "genus")
# 
# transformed_columns <- taxonomic_data %>%
#   mutate_at(vars(all_of(columns_to_process)), list(~ na.locf(.))) %>%
#   select(all_of(columns_to_process))
# 
# trait_data <- cbind(transformed_columns, taxonomic_data)


# trait_data <- condensed_traits %>%
#   select(all_of(columns_to_process)) %>%
#   t() %>%
#   na.locf() %>%
#   t() %>%
#   unselect()
