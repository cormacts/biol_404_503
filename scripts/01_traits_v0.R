#############
# BIOL 403 503 Project
# Script: 01_traits_v0
# By: Lydia, Elias
# 2024-03-18
# This script will organize functional traits based on the data from Weigel (2022) paper
#############

traits <- read.csv("data/raw/msystems.01422-21-s0005.csv")

## Condensing categories of traits into single columns

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

## Creating a new condensed dataframe 
# Removed MAG name, other traits not relevant to project

condensed_traits <- traits[1:66,c("phylum","class","order","family","genus","dissimilatory_nitrate_reduction",
                              "assimilatory_nitrate_reduction","denitrification","nitrogen_fixation","nitrification")]
