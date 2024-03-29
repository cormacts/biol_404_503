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

png(file = "figures/mag_taxaplot.png")
sp_stackplot
dev.off()

##### stats ####

#creat dataframes from each dataset (mag data and 16s data) with the proportions of nitrogen cyclers

s_proportions <- filter(plot_sp_proportions,species == "Nereocystis", nitrogen_cycling == "proportion_Y_22") 
  
relabun_Y <- filter(relabun_data, nitrogen_cycling == "YES")

proportion_16s <- s_proportions$proportions_22
proportion_mag <- relabun_data$relabun
# Perform Welch's t-test
ttest_mag_16s <- t.test(proportion_16s, proportion_mag, var.equal = FALSE)

# Print the result
print(ttest_mag_16s)


# Didn't check our assumptions but added this in case, also p<0.5
utest_mag_16s <- wilcox.test(proportion_16s, proportion_mag)
print(utest_mag_16s)

#t = 2.765, df = 40.312, p-value = 0.008546
# there is a statistically significant difference between the relative abundances of nitrogen cyclers in
# the metagenomic data vs the 16S data for community composition of functional traits. Therefore we can
# reject the null hypothesis that these two groups are the same. 

#This means that the functional trait composition
# data between these two groups cannot be compared because there are two many unknowns in the 16s data




