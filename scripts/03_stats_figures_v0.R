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
library(ggpubr)
library(rstatix)

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

## Checking if it worked
# View(species_data)
# summary(species_data)
# levels(species_data$nitrogen_cycling)

## Summing the asv_abundances based off nitrogen cycling for each species
sum_sp_data <- species_data %>%
  group_by(description, nitrogen_cycling) %>%
  summarise_at(vars(asv_abundance),
               list(sum_abundance = sum))
## double checking if the total asv_abundance values checks out:
##sum(species_data$asv_abundance) yes it did check out 


## Plan for Code:
## initial statistical analysis: 
## Individual t-test for between species and between location

# ## t-test: ## variable names subject to change later on
## getting the sums for nitrogen cycling microbes on each kelp species:
# mac_sum <- sum_sp_data$sum_abundance[sum_sp_data$nitrogen_cycling == "Yes" &
#                                      sum_sp_data$description == "Macrocystis"]
# 
# ner_sum <- sum_sp_data$sum_abundance[sum_sp_data$nitrogen_cycling == "Yes" &
#                                        sum_sp_data$description == "Nereocystis"]
## Okay, above method isn't really working for t-test purposes...


## Filtering the data to focus specifically on nitrogen cycling species
species_data%>%
  filter(nitrogen_cycling == "Yes") -> ndata

## Filtering data to focus on nitrogen cycling microbes for kelp sites
blade_data%>%
  filter(nitrogen_cycling == "Yes") -> bdata


## Shapiro test:
## Between kelp species
ndata %>%
  group_by(description) %>%
  shapiro_test(asv_abundance)
## Output from above, p-value < 0.05 which implies distribution of data is significantly
## different from normal distribution, cannot assume normality

## Between meristem and blade tip
bdata %>%
  group_by(sample_type) %>%
  shapiro_test(asv_abundance)
## Output from above, p-value < 0.05 which implies distribution of data is significantly
## different from normal distribution, cannot assume normality

## Levene's test for homogeneity of variances:
## For between kelp species:
sp_lt <- leveneTest(asv_abundance ~ description, ndata)
print(sp_lt)
## p-value is greater than 0.05, not enough evidence to reject null hypothesis

## For between meristem and blade tip:
b_lt <- leveneTest(asv_abundance ~ sample_type, bdata)
print(b_lt)
## p-value is less than 0.05, therefore we can reject the null hypothesis, can apply two-sample t-test


## Between species (sp)test: mean number of nitrogen fixing microbe species
## Using a Mann-Whitney U test as our data does not meet the assumptions of the t-test
sp_wctest <- wilcox.test(asv_abundance ~ description, ndata)
print(sp_wctest)
## p-value is below 0.05, therefore we reject the null hypothesis

## Between kelp location (klc) t-test: mean number of nitrogen fixing microbe species
b_ttest <- t.test(asv_abundance ~ sample_type, bdata)
print(b_ttest)
# p-value is below 0.05, therefore we reject the null hypothesis and consider the alternative hypothesis

## Plan for Figures:
## Bar plots: possible 2 graphs -> proportion
## Y - axis: relative abundance
## X- axis: nitrogen fixing (yes/no), unknown

## stack bar chart: (skeleton -> currently looking at ASV abundance)
## bar chart colours subject to change, I (Annie) just think these are cute and visible

## Stack bar graph looking at asv abundance of the microbes involved in nitrogen cycling,
## -not involved in nitrogen cycling and unknown. Comparison between abundances on Macrocystis and Nereocystis
sp_stackplot <- ggplot(species_data, aes(fill=nitrogen_cycling, y=asv_abundance, x=description)) + 
                    geom_bar(position="stack", stat="identity")+
                    labs(x = "Species", y = "ASV Abundance", color = "Nitrogen Cycling") +
                    theme(strip.text = element_text(face = "italic"),
                          axis.text.x = element_text(colour = "grey20", size = 12)) +
                    theme_bw() +
                    scale_fill_manual(values = c("lightblue", "yellow1", "violet"))
sp_stackplot
ggsave(file = "/figures/species_stackplot.PDF", plot = sp_stackplot, dpi = 500, units = "mm", width = 150, height = 100)

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


## Plans for species diversity plots: currently unsure how this will look at the moment
## Y-axis:
## X-axis:
## wrap:
