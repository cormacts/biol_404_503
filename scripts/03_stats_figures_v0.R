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
#install.packages("hrbrthemes")
library(hrbrthemes)
#install.packages("viridis")
library(viridis)

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

## Revising the dataframe used for the statistical tests:
new_sum_sp_data <- species_data %>%
  group_by(Row.names, nitrogen_cycling, description) %>%
  summarise_at(vars(asv_abundance),
               list(sum_abundance = sum))

# Pivoting so that the Y/N/Unknown abundances are columns with one row per sample
pivoted_sum_sp_data <- new_sum_sp_data %>%
  pivot_wider(names_from = nitrogen_cycling,values_from = sum_abundance)

### Adding columns calculating: 
# 1. Total abundance of all ASVs that are taxa covered by Weigel 2022 (functional traits paper)
pivoted_sum_sp_data$total_abundance_Weigel2022 = pivoted_sum_sp_data$Yes+pivoted_sum_sp_data$No
# 2. Total abundance of all ASVs covered by Weigel 2019 (original paper)
pivoted_sum_sp_data$total_abundance_Weigel2019 = pivoted_sum_sp_data$Yes+pivoted_sum_sp_data$No+pivoted_sum_sp_data$Unknown
# 3. Proportions for Y, N, and Unknown nitrogen cycling traits, out of either the 2022 or 2019 total
pivoted_sum_sp_data$proportion_Y_22 = pivoted_sum_sp_data$Yes/pivoted_sum_sp_data$total_abundance_Weigel2022
pivoted_sum_sp_data$proportion_N_22 = pivoted_sum_sp_data$No/pivoted_sum_sp_data$total_abundance_Weigel2022
pivoted_sum_sp_data$proportion_Y_19 = pivoted_sum_sp_data$Yes/pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_sp_data$proportion_N_19 = pivoted_sum_sp_data$No/pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_sp_data$proportion_NA_19 = pivoted_sum_sp_data$Unknown/pivoted_sum_sp_data$total_abundance_Weigel2019

## Filtering the data to focus specifically on nitrogen cycling species
## So we can compare directly the abundance of nitrogen cycling microbes
new_sum_sp_data%>%
  filter(nitrogen_cycling == "Yes") -> ndata

## Doing the same as above for blade data
new_sum_blade_data <- blade_data %>%
  group_by(Row.names, nitrogen_cycling,sample_type) %>%
  summarise_at(vars(asv_abundance),
               list(sum_abundance = sum))

# Pivoting so that the Y/N/Unknown abundances are columns with one row per sample
pivoted_sum_blade_data <- new_sum_blade_data %>%
  pivot_wider(names_from = nitrogen_cycling,values_from = sum_abundance)

# Adding columns calculating 
# 1. Total abundance of all ASVs that are taxa covered by Weigel 2022 (functional traits paper)
pivoted_sum_blade_data$total_abundance_Weigel2022 = pivoted_sum_blade_data$Yes+pivoted_sum_blade_data$No
# 2. Total abundance of all ASVs covered by Weigel 2019 (original paper)
pivoted_sum_blade_data$total_abundance_Weigel2019 = pivoted_sum_blade_data$Yes+pivoted_sum_blade_data$No+pivoted_sum_blade_data$Unknown
# 3. Proportions for Y, N, and Unknown nitrogen cycling traits, out of either the 2022 or 2019 total
pivoted_sum_blade_data$proportion_Y_22 = pivoted_sum_blade_data$Yes/pivoted_sum_blade_data$total_abundance_Weigel2022
pivoted_sum_blade_data$proportion_N_22 = pivoted_sum_blade_data$No/pivoted_sum_blade_data$total_abundance_Weigel2022
pivoted_sum_blade_data$proportion_Y_19 = pivoted_sum_blade_data$Yes/pivoted_sum_blade_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_N_19 = pivoted_sum_blade_data$No/pivoted_sum_blade_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_NA_19 = pivoted_sum_blade_data$Unknown/pivoted_sum_blade_data$total_abundance_Weigel2019

## Filtering data to focus on nitrogen cycling microbes for kelp sites
## So we can compare directly the abundance of nitrogen cycling microbes
new_sum_blade_data%>%
  filter(nitrogen_cycling == "Yes") -> bdata

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


## We want to do a two-sample comparison to look at possible differences in microbe functionality between kelp species
## In this case, we want to examine the differences between nitrogen cycling abundance found on the two kelp

## Testing for assumptions below --------------------
## Shapiro test:
## Between kelp species
ndata %>%
  group_by(description) %>%
  shapiro_test(sum_abundance)
## Output from above, p-value > 0.05 which implies distribution of data is under normal distribution
## Can assume normality

## Between meristem and blade tip
bdata %>%
  group_by(sample_type) %>%
  shapiro_test(sum_abundance)
## Output from above, p-value > 0.05 which implaies distribution of data is under normal disribution
## Can assume normality

## Levene's test for homogeneity of variances:
## For between kelp species:
sp_lt <- leveneTest(sum_abundance ~ description, ndata)
print(sp_lt)
## p-value is less than 0.05, therefore we can reject the null hypothesis: We can continue forward using a two-sample t-test

## For between meristem and blade tip:
b_lt <- leveneTest(sum_abundance ~ sample_type, bdata)
print(b_lt)
## p-value is greater than 0.05, therefore we cannot reject the null hypothesis:
## will perform a Mann-Whitney U test

## --------------------------------------------------

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
ggsave(file = "figures/species_stackplot.PDF", plot = sp_stackplot, dpi = 500, units = "mm", width = 150, height = 100)

# Observations from the plot: 
# Looks as though Macrocystis has a larger abundance of nitrogen cycling microbes, which contradicts what we initially hypothesied, 
# however it is possible that when comparing proportion of nitrogen cycling microbes
# to the rest of the population, we may find different results. 
# In addition, it is important to keep in mind that a large abundance of microbes 
# found in both species are currently classified as "Unknown".
# So there may be more nitrogen cycling species within that group. 

## Stack bar graph looking at asv abundance of the microbes involved in nitrogen cycling,
## -not involved in nitrogen cycling and unknown. Comparison between abundances on meristem and blade tip
blade_stackplot <- ggplot(blade_data, aes(fill=nitrogen_cycling, y=asv_abundance, x=sample_type)) + 
                      geom_bar(position="stack", 
                               stat="identity") +
                      labs(x = "Sample Type", y = "ASV Abundance", color = "Nitrogen Cycling") +
                      theme(strip.text = element_text(face = "italic"),
                            axis.text.x = element_text(colour = "grey20", size = 12)) +
                      theme_bw() +
                      scale_fill_manual(values = c("lightblue", "yellow1", "violet"))

blade_stackplot
ggsave(file = "figures/blade_stackplot.PDF", plot = blade_stackplot, dpi = 500, units = "mm", width = 150, height = 100)

## Observations of plot:
## There is a vast difference between abundance of nitrogen cycling on the kelp blade tip and meristem.
## In contrast to our hypothesis, there is a higher abundance of nitrogen cycling microbes on the blade tip than the kelp meristem.
## There are also in general, way more microbes on the blade tip than the meristem. Again keeping in mind that there are still many
## microbes classified as "Unknown" and these could include nitrogen cycling microbes.

## Plans for box-plots:
## Y-axis: Either Proportion (in decimals) or Percentage (%)
## X-axis: Group
## Plots: Proportions for each individual sample (box plot will be the average?)

## Box plot comparing proportions of nitrogen cycling microbes between Macrocystis and Nereocystis
plot_blade_proportions %>%
  filter(nitrogen_cycling == "proportion_Y_22") -> b_plot_data




## Box plot comparing proportions of nitrogen cycling microbes between kelp meristem and blade tip samples
## Filtering for nitrogen cycling microbes
plot_blade_proportions %>%
  filter(nitrogen_cycling == "proportion_Y_22") -> b_plot_data

## Code to plot the graph
b_plot_data %>%
  ggplot( aes(x=blade_location, y=proportions_22, fill=nitrogen_cycling)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Boxplot of Nitrogen Cycling Proportions") +
  xlab("")
