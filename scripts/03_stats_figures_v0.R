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


cleanwrp = read.csv("data/processed/cleanwrp_dataframe.csv")
View(cleanwrp)

## Testable hypotheses:
# The functional roles of microbiota associated with Macrocystis and Nereocystis differ due to the distinct life histories of these two kelp species:
# Genes associated with nitrogen reduction will be more prevalent in Nereocystis due to its high rate of growth, being an annual kelp species.

# The functional roles of microbiota associated with younger, meristematic blade tissue differ from those associated with older, apical blade tissue in Nereocystis luetkeana:
# Genes associated with nitrogen reduction will be present at a higher proportion in microbial communities at the meristem as this is where growth is occurring in the kelp individual, creating biological available nitrogen in these regions.

## Plan for Code:
## initial statistical analysis: 
## Individual t-test for between species and between locations

## checking assumptions for individual t-test
## Levene's test to test for equal variance


## t-test: ## variable names subject to change later on
## between species (sp) t-test: mean number of nitrogen fixing microbe species
sp_t_test <- t.test(nereocystis_data, macrocystis_data)
print(spp_t_test)

## between kelp location (klc) t-test: mean number of nitrogen fixing microbe species
klc_t_test <- t.test(meristem_data, blade_data)
print(klc_t_test)


## Plan for Figures:
## Bar plots:
## variable subject to change:
p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")

p

