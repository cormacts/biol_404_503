#######
# BIOL 403 - 503 Project
# Figure Generation & Statistical Analysis
# Annie Wang

# Tentative methodology for analysis & data visualization:
# Importing required packages
# Reading in data files (filtered and rarefied)

library(tidyverse)
library(car)
library(phyloseq)
library(ggplot2)
library(forcats)
library(dplyr)

## Using dataframe made via previous scripts

## Testable hypotheses:
# The functional roles of microbiota associated with Macrocystis and Nereocystis differ due to the distinct life histories of these two kelp species:
# Genes associated with nitrogen reduction will be more prevalent in Nereocystis due to its high rate of growth, being an annual kelp species.

# The functional roles of microbiota associated with younger, meristematic blade tissue differ from those associated with older, apical blade tissue in Nereocystis luetkeana:
# Genes associated with nitrogen reduction will be present at a higher proportion in microbial communities at the meristem as this is where growth is occurring in the kelp individual, creating biological available nitrogen in these regions.

## Making Y/N/NA into more comprehensive categorical variables
## all missing data (N/A) values will be assigned "Unknown"
species_data <- species_data %>% 
  mutate(
    nitrogen_cycling = recode_factor(nitrogen_cycling,
                          "N" = "No",
                          "Y" = "Yes"),
    nitrogen_cycling = fct_explicit_na(nitrogen_cycling, "Unknown"))

blade_data <- blade_data %>% 
  mutate(
    nitrogen_cycling = recode_factor(nitrogen_cycling,
                                     "N" = "No",
                                     "Y" = "Yes"),
    nitrogen_cycling = fct_explicit_na(nitrogen_cycling, "Unknown"))

View(species_data)
summary(species_data)
levels(species_data$nitrogen_cycling)

mac_data <- species_data %>%
  group_by(nitrogen_cycling, description) %>%
  summarise(mac_sum = sum(asv_abundance))



View(blade_data)

## Aggregating data
## Summing of asv_abundance for nitrogen cyclers, non-N cyclers, and unknowns
## using a test file so i don't mess up the actual dataframe 

sp_test <- sp_test %>%
  group_by(nitrogen_cycling) %>%
  summarise(abundance = sum(asv_abundance))

# species_data %>%
#   filter(description == "Macrocystis") -> mac_test
# 
# species_data %>%
#   filter(description == "Nereocystis") -> ner_test
# 
# summary(ner_test)
# summary(mac_test)
#   
# mac_prop <- mac_test %>%
#   group_by(nitrogen_cycling) %>%
#   summarise(abundance = sum(asv_abundance))
#   
# summary(mac_prop)


#ner_test %>%
  group_by(location) %>%
  summarise(abundance = sum(asv_abundance))


View(mac_prop)

## Plan for Code:
## initial statistical analysis: 


## STATS BELOW ARE SUBJECT TO CHANGE, THIS IS JUST ME (ANNIE) BEING LOST 
## AND WORKING THROUGH SOME THINGS
## anova test:
## between species that are nitrogen fixing, not nitrogen fixing, and unknowns.
# a_sp = aov(asv_abundance ~ nitrogen_cycling, data = species_data)
# summary(a_sp)
# plot(a_sp)
# ## anova conditions are not met
# ##Kruskal-Wallis Test
# kwt_sp = kruskal.test(asv_abundance ~ nitrogen_cycling, data = species_data)
# summary(kwt_sp)
# ## Pairwise wilcoxon
# pairwise.wilcox.test(species_data$asv_abundance, species_data$nitrogen_cycling, p.adjust.method = "BH")


## Individual t-test for between species and between location

## checking assumptions for individual t-test
## Levene's test to test for equal variance


# ## t-test: ## variable names subject to change later on
# ## between species (sp) t-test: mean number of nitrogen fixing microbe species
sp_t_test <- t.test( macrocystis_data)
print(spp_t_test)

## between kelp location (klc) t-test: mean number of nitrogen fixing microbe species
klc_t_test <- t.test(meristem_data, blade_data)
print(klc_t_test)


## Plan for Figures:
## Bar plots: possible 2 graphs -> proportion
## Y - axis: relative abundance
## X- axis: nitrogen fixing (yes/no), unknown

## stack bar chart: (skeleton -> currently looking at ASV abundance)
## bar chart colours subject to change, I (annie) just think these are cute and visible

## Stack bar graph looking at asv abundance of the microbes involved in nitrogen cycling,
## -not involved in nitrogen cycling and unknown. Comparison between abundances on Macrocystis and Nereocystis
ggplot(species_data, aes(fill=nitrogen_cycling, y=asv_abundance, x=description)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Species", y = "ASV Abundance", color = "Nitrogen Cycling") +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(colour = "grey20", size = 12)) +
  theme_bw() +
  scale_fill_manual(values = c("lightblue", "yellow1", "violet"))

## Stack bar graph looking at asv abundance of the microbes involved in nitrogen cycling,
## -not involved in nitrogen cycling and unknown. Comparison between abundances on meristem and blade tip
ggplot(blade_data, aes(fill=nitrogen_cycling, y=asv_abundance, x=sample_type)) + 
  geom_bar(position="stack", 
           stat="identity") +
  labs(x = "Sample Type", y = "ASV Abundance", color = "Nitrogen Cycling") +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(colour = "grey20", size = 12)) +
  theme_bw() +
  scale_fill_manual(values = c("lightblue", "yellow1", "violet"))

