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
mag_abundance <- read.csv("data/raw/metag_abun.csv")
# lowest taxonomic unit column
mag_data_na <- mag_data_df %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(Genus), Genus,
                               if_else(!is.na(Family), Family,
                                       if_else(!is.na(Order), Order,
                                               if_else(!is.na(Class), Class,
                                                       if_else(!is.na(Phylum), Phylum, NA_character_))))))
#replace empty with na function
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}
#use function
mag_data_df <- mag_data %>%
  mutate(across(everything(), replace_empty_with_na))

#merge dataframes
mag_merged <- left_join(mag_data_na, trait_data, by = "lowest_rank")

#select columns for condensed data frame
mag_data_all <- mag_sum %>%
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
mag_total <- pivot_mag %>% group_by(Sample) %>% mutate(sample_total = sum(count))
# create cycling abundance column grouping by sample and cycling and calculate relative abundance
mag_total_relabun <- mag_total %>%  group_by(Sample, nitrogen_cycling) %>% mutate(nitro_total = sum(count)) %>% mutate(relabun = nitro_total/total)

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

png(file = "figures/mag_proportion_taxaplot.png")
sp_stackplot
dev.off()
