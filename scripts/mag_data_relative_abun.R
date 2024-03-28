#By: Cormac 
#calculating and plotting relative abundance of different microbes
#for nitrogen_cycling within different samples
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(vegan)
library(plyr)
library(qualpalr)
library(ggpubr)

mag_abundance <- read.csv("data/raw/msystems.01422-21-s0003.csv") 
trait_data <- read.csv("trait_data_w2022.csv")


#replace empty with na function
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}
#use function
mag_data_na <- mag_abundance %>%
  mutate(across(everything(), replace_empty_with_na))
# lowest taxonomic unit column
mag_data_df <- mag_data_na %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(Genus), Genus,
                               if_else(!is.na(Family), Family,
                                       if_else(!is.na(Order), Order,
                                               if_else(!is.na(Class), Class,
                                                       if_else(!is.na(Phylum), Phylum, NA_character_))))))

#merge dataframes
mag_merged <- left_join(mag_data_df, trait_data, by = "lowest_rank")

mag_merged <- mag_merged %>% 
  mutate(
    nitrogen_cycling = recode_factor(nitrogen_cycling,
                                     "N" = "No",
                                     "Y" = "Yes"),
    nitrogen_cycling = fct_explicit_na(nitrogen_cycling, "Unknown"))
#select columns for condensed data frame
mag_data_all <- mag_merged %>%
  select(
    MAG_name, lowest_rank,
   B1_Squaxin_2019:B7_Tatoosh_2018,
    nitrogen_cycling)
#make dataframe longer and put sites/samples in a column, create a count column
pivot_mag <- mag_data_all %>% 
  pivot_longer(
    cols = B1_Squaxin_2019:B7_Tatoosh_2018, 
    names_to = "Sample", 
    values_to = "count")


# create total abundance column grouping by sample
mag_total <- pivot_mag %>% group_by(Sample) %>% mutate(total = sum(count))
# create cycling abundance column grouping by sample and cycling and calculate relative abundance
mag_total_relabun <- mag_total %>%  
  group_by(Sample, nitrogen_cycling) %>% 
  mutate(nitro_total = sum(count)) %>% mutate(relabun = nitro_total/total)

#select for columns, pick out unique rows 
relabun_data <- mag_total_relabun %>%
  select(nitrogen_cycling, Sample, relabun) %>% 
  unique()

#plot
sp_stackplot <- ggplot(relabun_data, aes(fill=nitrogen_cycling, y=relabun, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Samples", y = "Relative Abundance", color = "Nitrogen Cycling") +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(colour = "grey20", size = 12)) +
  theme_bw() +
  scale_fill_manual(values = c("lightblue", "yellow1", "violet"))
sp_stackplot

# stats

plot_sp_proportions

## Levene's test for homogeneity of variances:
## For between mag data (2022) and 16s data (2019):
sp_lt <- leveneTest(sum_abundance ~ description, ndata)
print(sp_lt)
## p-value > 0.05, therefore we cannot reject the null hypothesis: 
## Cannot use a two-sample t-test
## Will perform a Mann-Whitney U test

## For between meristem and blade tip:
b_lt <- leveneTest(sum_abundance ~ sample_type, bdata)
print(b_lt)
## p-value is greater than 0.05, therefore we cannot reject the null hypothesis:
## Cannot use a two-sample t-test
## Will perform a Mann-Whitney U test
