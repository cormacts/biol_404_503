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

## Using dataframe made via previous scripts

## Testable hypotheses:
# The functional roles of microbiota associated with Macrocystis and Nereocystis differ due to the distinct life histories of these two kelp species:
# Genes associated with nitrogen reduction will be more prevalent in Nereocystis due to its high rate of growth, being an annual kelp species.

# The functional roles of microbiota associated with younger, meristematic blade tissue differ from those associated with older, apical blade tissue in Nereocystis luetkeana:
# Genes associated with nitrogen reduction will be present at a higher proportion in microbial communities at the meristem as this is where growth is occurring in the kelp individual, creating biological available nitrogen in these regions.

## Plan for Code:
## initial statistical analysis: 
## Individual t-test for between species and between locations

## Making Y/N/NA into more comprehensive categorical variables
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
View(blade_data)

## checking assumptions for individual t-test
## Levene's test to test for equal variance


## STATS BELOW ARE SUBJECT TO CHANGE, THIS IS JUST ME (ANNIE) BEING LOST 
## AND WORKING THROUGH SOME THINGS
## anova test:
## between species that are nitrogen fixing, not nitrogen fixing, and unknowns.
a_sp = aov(asv_abundance ~ nitrogen_cycling, data = species_data)
summary(a_sp)
plot(a_sp)

## anova conditions are not met
##Kruskal-Wallis Test
kwt_sp = kruskal.test(asv_abundance ~ nitrogen_cycling, data = species_data)
summary(kwt_sp)

## Pairwise wilcoxon
pairwise.wilcox.test(species_data$asv_abundance, species_data$nitrogen_cycling, p.adjust.method = "BH")

## t-test: ## variable names subject to change later on
## between species (sp) t-test: mean number of nitrogen fixing microbe species
sp_t_test <- t.test(nereocystis_data, macrocystis_data)
print(spp_t_test)

## between kelp location (klc) t-test: mean number of nitrogen fixing microbe species
klc_t_test <- t.test(meristem_data, blade_data)
print(klc_t_test)


## Plan for Figures:
## Bar plots: possible 2 graphs -> proportion
## Y - axis: relative abundance
## X- axis: nitrogen fixing (yes/no), unknown

## stack bar chart: 

## variable subject to change:
plot(asv_abundance ~ nitrogen_cycling, data = species_data)

## general outline: still need to sum species abundance
ggplot(data = species_data, aes(x = nitrogen_cycling, y = asv_abundance, fill = description)) +
  geom_bar(stat = "summary",
           fun.data = "mean_se",
           position = "dodge") +
  geom_errorbar(stat = "summary",
                fun.data = "mean_se",
                position = position_dodge(width = 0.9),
                width = 0.2) +
  labs(x = "", 
       y = " ") +
  theme(strip.text = element_text(face = "italic"))
