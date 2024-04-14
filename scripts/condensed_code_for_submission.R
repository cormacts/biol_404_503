### BIOL 403/503 project ###
# Condensed code for submission

### Initial setup ###

# Loading Libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(plyr)
library(qualpalr)
library(ggpubr)
library(dplyr)
library(car)
library(ggplot2)
library(forcats)
library(rstatix)
#install.packages("hrbrthemes")
library(hrbrthemes)
#install.packages("viridis")
library(viridis)

# Reading in data from Weigel 2019
wrp = readRDS("data/raw/W2019_and_RP2022_unfiltered_phyloseq.RDS")

### Filtering ###
# Code adapted from BIOL 403/503 Lab 2

## 1. Removing off target taxa (we are considering only Bacteria and Archaea)
wrp = subset_taxa(
  wrp,
  domain != "Unassigned" &
    order != "Chloroplast" &
    order != "Mitochondria" &
    family != "Mitochondria" &
    domain != "Eukaryota"
)

## 2. Removing samples with low number of reads
sampleSums(wrp) # prints sample sums
# Adding the sample sums of the reads to meta data column
wrp@sam_data$sample_sums_unfiltered = as.numeric(sample_sums(wrp))
# We see a relatively continuous increase in read numbers per sample, so will use 800 as the cutoff.
# Removing all samples with less than 800 reads
wrp.high <- prune_samples(sample_sums(wrp) >= 800, wrp)
# Making a new object with only the samples  lost (less than 800 reads)
wrp.below800 <- prune_samples(sample_sums(wrp) < 800, wrp)
# Getting the metadata out of phyloseq for low reads object
wrp.below800 = as.matrix(wrp.below800@sam_data)
# Writing file to know which samples were lost
write.csv(wrp.below800,
          "data/processed/wrp_samples_less_than_800.csv")

## 3. Removing ASVs that are low frequency
# Extracting otu dataframe (asv table) from phyloseq object
otutab <- as.data.frame(t(as.matrix(otu_table(wrp.high@otu_table))))
# Calculating the sum of each row in the otutab
otutab$asv_abundance = rowSums(otutab)
# Finding the minimum value of asv_abundance
min(otutab$asv_abundance)
# Removing ASVs with low frequency using 100 as threshold
otu.pruned = subset(otutab, otutab$asv_abundance >= 100)
# Confirming the new minimum ASV abundance value is at the threshold
min(otu.pruned$asv_abundance)
# Removing ASV_abundance column
widthotu = ncol(otu.pruned) # finding width
otu.pruned = otu.pruned[, -c(widthotu)] # keep everything except the last columns

## 4. Removing infrequent ASVs over samples
# Creating a function that counts the number of occurrence along rows
ASVoccur = function(x) {
  return(sum(x > 0))
}
# Calculating the occurrence of each ASV in the dataframe
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)
# Investigating what the ASV occurance counts look like
summary(otu.pruned$asv_occur_count)
# Removing ASVs found two or less times
otu.highfreq = subset(otu.pruned, otu.pruned$asv_occur_count > 2)
# Confirming filtering worked
summary(otu.highfreq$asv_occur_count)
# Removing the asv_occur_count column
otu.highfreq = otu.highfreq[, -c(widthotu)]

## 5. Data de-noising
otu.clean <- mutate_all(otu.highfreq, funs(ifelse(. < 3, 0, .)))

### Making a new phyloseq object of cleaned data ###

cleanwrp = phyloseq(
  sample_data(wrp.high@sam_data),
  tax_table(wrp.high@tax_table),
  otu_table(as.matrix(otu.clean), taxa_are_rows = TRUE)
)


## Creating the dephyloseq function: breaks down the phyloseq object into more manageable data
dephyloseq = function(phylo_obj) {
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  ## how many metadata columns you have
  metacols = ncol(meta) + 1
  ## get out the otu table
  ## if your metadta is empty after running this, you need to use
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by = 0)
  ## get out the taxonomy file
  tax = as.data.frame(phylo_obj@tax_table)
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var = "asv_name")
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table
  mo = mo %>% pivot_longer(
    cols = -c(1:metacols),
    names_to = "asv_name",
    values_to = "asv_abundance"
  )
  ## Join the metadata and otu table with the taoxnomy table
  mot = full_join(mo, tax)
  ## Specify the output for the dephyloseq function
  output = mot
}

## We are not rarefying our data at this stage, since we will work with proportions later.

## Getting cleanwrp data out of phyloseq and into a dataframe
community_data = dephyloseq(cleanwrp)

## Writing to a file
# write.csv(community_data, "data/processed/community_dataframe.csv")


### Setting up the functional trait dataset ###

## Importing Data
traits <- read.csv("data/raw/msystems.01422-21-s0005.csv")

## Condensing all nitrogen cycling-related traits into one column
traits <- traits %>%
  mutate_at(vars(starts_with(
    c(
      "Dissimilatory_nitrate_reduction",
      "Nitrogen_fixation",
      "Nitrification"
    )
  )),  ~ coalesce(., "NA")) %>%
  mutate(nitrogen_cycling = if_else(rowSums(select(., starts_with(
    c(
      "Dissimilatory_nitrate_reduction",
      "Nitrogen_fixation",
      "Nitrification"
    )
  )) == "Y", na.rm = TRUE) > 0, "Y", ""))

## Creating a new dataframe with only relevant variables
condensed_traits <-
  traits[1:66, c("MAG_name",
                 "phylum",
                 "class",
                 "order",
                 "family",
                 "genus",
                 "nitrogen_cycling")]

## Finding lowest taxonomic data to apply to other dataset
# Creating function to remove blank strings and replace with NAs
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}
# Applying the above function to our dataset
condensed_traits <- condensed_traits %>%
  mutate(across(everything(), replace_empty_with_na))
# Creating a new column with the lowest taxonomic rank available
trait_data <- condensed_traits %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(genus), genus,
                               if_else(
                                 !is.na(family), family,
                                 if_else(!is.na(order), order,
                                         if_else(
                                           !is.na(class), class,
                                           if_else(!is.na(phylum), phylum, NA_character_)
                                         ))
                               )))
# The previous function will have placed NAs in the nitrogen_cycling column where there were empty spaces.
# Creating a function to replace those NAs with "N" (= no nitrogen cycling traits)
replace_na_with_n <- function(x) {
  x[is.na(x)] <- "N"
  return(x)
}
# Applying function to our dataset
trait_data <- trait_data %>%
  mutate(nitrogen_cycling = replace_na_with_n(nitrogen_cycling))

# To deal with repeat cases of lowest_rank, we are choosing to merge all observations of a given taxon into a single row.
# This will put a "Y" for nitrogen_cycling if any of the observations in that group had a "Y".
trait_data <- trait_data %>% arrange(desc(nitrogen_cycling))
trait_data <- trait_data %>%
  group_by(lowest_rank) %>%
  slice(1) %>%
  ungroup()


### Merging the two datasets ###

# community_data contains data relating to community composition and shifts based on treatment (Weigel et al., 2019)
# trait_data contains microbial taxa and their functional traits (Weigel et al., 2022)

## We will join these two dataframes using left_join(), to keep all observations in the left (first) dataframe.
all_merged <-
  left_join(community_data, trait_data, by = join_by("species" == "lowest_rank"))
# Renaming the genus column to reflect that the lowest available taxonomic rank is now stored there
names(all_merged)[names(all_merged) == "species"] <- "lowest_rank"

## Creating two sets of filtered data, one for each set of conditions we are planning to compare (meristem versus tip / Macrocystis vs. Nereocystis host)
# Either question can only be asked using data from specific sites - need to filter for these locations.

# In species_data - only want Macrocystis and Nereocystis, and locations with both species present
both_species_locations <-
  c("Bullman",
    "Cape_Johnson",
    "Destruction_Island",
    "Koitlah",
    "Sekiu")
species_data <- all_merged %>%
  select(
    Row.names,
    description,
    host_common_name,
    location,
    sample_type_user,
    asv_abundance,
    lowest_rank,
    nitrogen_cycling
  ) %>%
  filter(description == "Macrocystis" |
           description == "Nereocystis") %>%
  filter(location %in% both_species_locations)

# In blade_data - only want Nereocystis base and Nereocystis tip (whats the difference between "tip" and "tip_swab")
base_and_tip <-
  c(
    "Nereocystis_tip",
    "Nereocystis_base",
    "Nereocystis_tip_swab",
    "Nereocystis_base_swab"
  )
blade_data <- all_merged %>%
  select(
    Row.names,
    description,
    host_common_name,
    location,
    sample_type_user,
    asv_abundance,
    lowest_rank,
    nitrogen_cycling
  ) %>%
  filter(description %in% base_and_tip) %>%
  filter(location == "Tatoosh" |
           location == "Tatoosh_Main_Beach") %>%
  mutate(sample_type = if_else(
    grepl("tip", description),
    "tip",
    if_else(grepl("base", description), "meristem", NA_character_)
  ))

## Making Y/N/NA into more comprehensive categorical variables
## all missing data (N/A) values will be assigned "Unknown"
species_data <- species_data %>%
  mutate(
    nitrogen_cycling = recode_factor(nitrogen_cycling,
                                     "N" = "No",
                                     "Y" = "Yes"),
    nitrogen_cycling = fct_explicit_na(nitrogen_cycling, "Unknown")
  )

blade_data <- blade_data %>%
  mutate(
    nitrogen_cycling = recode_factor(nitrogen_cycling,
                                     "N" = "No",
                                     "Y" = "Yes"),
    nitrogen_cycling = fct_explicit_na(nitrogen_cycling, "Unknown")
  )


### Calculating proportions of microbes present with nitrogen cycling-related traits ###

## Grouping by sample (based on the IDs stored in Row.names)
sum_sp_data <- species_data %>%
  group_by(Row.names, nitrogen_cycling, description) %>%
  summarise_at(vars(asv_abundance),
               list(sum_abundance = sum))
sum_blade_data <- blade_data %>%
  group_by(Row.names, nitrogen_cycling, sample_type) %>%
  summarise_at(vars(asv_abundance),
               list(sum_abundance = sum))

## Pivoting so that the Y/N/Unknown abundances are columns with one row per sample
pivoted_sum_sp_data <- sum_sp_data %>%
  pivot_wider(names_from = nitrogen_cycling, values_from = sum_abundance)
pivoted_sum_blade_data <- sum_blade_data %>%
  pivot_wider(names_from = nitrogen_cycling, values_from = sum_abundance)

## Adding columns calculating:
# 1. Total abundance of all ASVs that are taxa covered by Weigel 2022 (functional traits paper)
pivoted_sum_sp_data$total_abundance_Weigel2022 = pivoted_sum_sp_data$Yes +
  pivoted_sum_sp_data$No
pivoted_sum_blade_data$total_abundance_Weigel2022 = pivoted_sum_blade_data$Yes +
  pivoted_sum_blade_data$No
# 2. Total abundance of all ASVs covered by Weigel 2019 (original paper)
pivoted_sum_sp_data$total_abundance_Weigel2019 = pivoted_sum_sp_data$Yes +
  pivoted_sum_sp_data$No + pivoted_sum_sp_data$Unknown
pivoted_sum_blade_data$total_abundance_Weigel2019 = pivoted_sum_blade_data$Yes +
  pivoted_sum_blade_data$No + pivoted_sum_blade_data$Unknown
# 3. Proportions for Y, N, and Unknown nitrogen cycling traits, out of either the 2022 or 2019 total
pivoted_sum_sp_data$proportion_Y_22 = pivoted_sum_sp_data$Yes / pivoted_sum_sp_data$total_abundance_Weigel2022
pivoted_sum_blade_data$proportion_Y_22 = pivoted_sum_blade_data$Yes / pivoted_sum_blade_data$total_abundance_Weigel2022
pivoted_sum_sp_data$proportion_N_22 = pivoted_sum_sp_data$No / pivoted_sum_sp_data$total_abundance_Weigel2022
pivoted_sum_blade_data$proportion_N_22 = pivoted_sum_blade_data$No / pivoted_sum_blade_data$total_abundance_Weigel2022
pivoted_sum_sp_data$proportion_Y_19 = pivoted_sum_sp_data$Yes / pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_Y_19 = pivoted_sum_blade_data$Yes / pivoted_sum_blade_data$total_abundance_Weigel2019
pivoted_sum_sp_data$proportion_N_19 = pivoted_sum_sp_data$No / pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_N_19 = pivoted_sum_blade_data$No / pivoted_sum_blade_data$total_abundance_Weigel2019
pivoted_sum_sp_data$proportion_NA_19 = pivoted_sum_sp_data$Unknown / pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_NA_19 = pivoted_sum_blade_data$Unknown /
  pivoted_sum_blade_data$total_abundance_Weigel2019

### Visualizing results ###

## Making taxaplot-like bar charts to show proportions

## Making dataframes in a format easy to plot
# Species comparison
pivoted_sum_sp_data$sample_ID = paste(pivoted_sum_sp_data$description,
                                      pivoted_sum_sp_data$Row.names)
plot_sp_proportions <-
  pivoted_sum_sp_data[, c("sample_ID", "proportion_Y_22", "proportion_N_22")]
plot_sp_proportions <- plot_sp_proportions %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling", values_to = "proportions_22")
plot_sp_proportions <- plot_sp_proportions %>%
  separate(sample_ID,
           c("species", "sample"),
           sep = " ",
           remove = TRUE)
# Blade location comparison
pivoted_sum_blade_data$sample_ID = paste(pivoted_sum_blade_data$sample_type,
                                         pivoted_sum_blade_data$Row.names)
plot_blade_proportions <-
  pivoted_sum_blade_data[, c("sample_ID", "proportion_Y_22", "proportion_N_22")]
plot_blade_proportions <- plot_blade_proportions %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling", values_to = "proportions_22")
plot_blade_proportions <- plot_blade_proportions %>%
  separate(sample_ID,
           c("blade_location", "sample"),
           sep = " ",
           remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
# Species comparison:
ggplot(plot_sp_proportions,
       aes(x = sample, y = proportions_22,
           fill = nitrogen_cycling)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("lightblue", "yellow1")) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(. ~ species, scales = "free", space = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(color = "black", size = 10),
    legend.text = element_text(size = 6),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank()
  ) +
  labs(y = "Relative abundance", x = "Sample", fill = "Nitrogen cycling")
# Blade location comparison:
ggplot(plot_blade_proportions,
       aes(x = sample, y = proportions_22,
           fill = nitrogen_cycling)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("lightblue", "yellow1")) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(. ~ blade_location, scales = "free", space = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(color = "black", size = 10),
    legend.text = element_text(size = 6),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank()
  ) +
  labs(y = "Relative abundance", x = "Sample", fill = "Nitrogen cycling")

## Repeating above procedure, but including all taxa covered by Weigel 2019

## Making dataframe to use for plotting
# Species comparison:
plot_sp_proportions_2 <-
  pivoted_sum_sp_data[, c("sample_ID",
                          "proportion_Y_19",
                          "proportion_N_19",
                          "proportion_NA_19")]
plot_sp_proportions_2 <- plot_sp_proportions_2 %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling", values_to = "proportions_19")
plot_sp_proportions_2 <- plot_sp_proportions_2 %>%
  separate(sample_ID,
           c("species", "sample"),
           sep = " ",
           remove = TRUE)
# Blade location comparison:
plot_blade_proportions_2 <-
  pivoted_sum_blade_data[, c("sample_ID",
                             "proportion_Y_19",
                             "proportion_N_19",
                             "proportion_NA_19")]
plot_blade_proportions_2 <- plot_blade_proportions_2 %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling", values_to = "proportions_19")
plot_blade_proportions_2 <- plot_blade_proportions_2 %>%
  separate(sample_ID,
           c("blade_location", "sample"),
           sep = " ",
           remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
# Species comparison:
ggplot(plot_sp_proportions_2,
       aes(x = sample, y = proportions_19,
           fill = nitrogen_cycling)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("lightblue", "violet", "yellow1")) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(. ~ species, scales = "free", space = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(color = "black", size = 10),
    legend.text = element_text(size = 6),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank()
  ) +
  labs(y = "Relative abundance", x = "Sample", fill = "Nitrogen cycling")
# Blade location comparison:
ggplot(plot_blade_proportions_2,
       aes(x = sample, y = proportions_19,
           fill = nitrogen_cycling)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("lightblue", "violet", "yellow1")) +
  guides(fill = guide_legend(ncol = 2)) +
  facet_grid(. ~ blade_location, scales = "free", space = "free") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(color = "black", size = 10),
    legend.text = element_text(size = 6),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank()
  ) +
  labs(y = "Relative abundance", x = "Sample", fill = "Nitrogen cycling")

## Making boxplot to compare proportions

## Filtering for nitrogen cycling microbe proportion for stat tests and boxplots:
sp_stat_data <- plot_sp_proportions %>%
  filter(nitrogen_cycling == "proportion_Y_22")
blade_stat_data <- plot_blade_proportions %>%
  filter(nitrogen_cycling == "proportion_Y_22")

## Box plot comparing proportions of nitrogen cycling microbes between Macrocystis and Nereocystis
## Code to plot graph
sp_boxplot <-
  ggplot(sp_stat_data, aes(x = species, y = proportions_22, fill = species)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",
              size = 0.4,
              alpha = 0.9) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title.x = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 12
    ),
    axis.title.y = element_text(hjust = 0.5, size = 12)
  ) +
  ggtitle("Boxplot of Nitrogen Cycling Proportions") +
  labs(x = "Kelp Species", y = "Proportion of Microbe Species")

png(file = "figures/sp_boxplot.png")
sp_boxplot
dev.off()

## Box plot comparing proportions of nitrogen cycling microbes between kelp meristem and blade tip samples

## Code to plot the graph
blade_boxplot <-
  ggplot(blade_stat_data,
         aes(x = blade_location, y = proportions_22, fill = blade_location)) +
  geom_boxplot() +
  labs(x = "Sample Location", y = "Proportion of Microbe Species") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",
              size = 0.4,
              alpha = 0.9) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11),
    axis.title.x = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 12
    ),
    axis.title.y = element_text(hjust = 0.5, size = 12)
  ) +
  ggtitle("Boxplot of Nitrogen Cycling Proportions") +
  labs(x = "Kelp Blade Locations", y = "Proportion of Microbe Species")

png(file = "figures/blade_boxplot.png")
blade_boxplot
dev.off()

### Statistical tests ###

## We want to do a two-sample comparison to look at possible differences in microbe functionality between kelp species and between blade location

## Testing for assumptions of a two-sample t-test:

## Shapiro test
## Between kelp species:
sp_stat_data %>%
  group_by(species) %>%
  shapiro_test(proportions_22)
# Output from above, p-value > 0.05 which implies distribution of data is under normal distribution
# Can assume normality

## Between meristem and blade tip:
blade_stat_data %>%
  group_by(blade_location) %>%
  shapiro_test(proportions_22)
# Output from above, p-value < 0.05, implies data is not under normal distribution
# Cannot assume normality as the two groups differ

## Levene's test for homogeneity of variances
## Between kelp species:
sp_lt <- leveneTest(proportions_22 ~ species, sp_stat_data)
print(sp_lt)
## p-value < 0.05, therefore we can reject the null hypothesis: There is not equal variance
## Cannot use a two-sample t-test --> will perform a Mann-Whitney U test

## For between meristem and blade tip:
b_lt <- leveneTest(proportions_22 ~ blade_location, blade_stat_data)
print(b_lt)
## p-value is greater than 0.05, therefore we cannot reject the null hypothesis
## Can use a two-sample t-test

## --------------------------------------------------

## Between species test: mean proportion of nitrogen-fixing microbe species
## Using a Mann-Whitney U test as our data does not meet the assumptions of the t-test
sp_wctest <- wilcox.test(proportions_22 ~ species, sp_stat_data)
print(sp_wctest)
## p-value > 0.05, therefore we cannot reject the null hypothesis

## Between kelp location (meristem and blade tip) test: mean proportion of nitrogen-fixing microbe species
## Using a Mann-Whitney U test as our data does not meet the assumptions of the t-test
b_wctest <-
  wilcox.test(proportions_22 ~ blade_location, blade_stat_data)
print(b_wctest)
## p-value < 0.05, therefore we can reject the null hypothesis and embrace the alternative hypothesis
