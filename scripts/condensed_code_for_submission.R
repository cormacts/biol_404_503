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

# Reading in data from Weigel 2019
wrp = readRDS("data/raw/W2019_and_RP2022_unfiltered_phyloseq.RDS")

# Filtering

# 1. Removing off target taxa (we are considering only Bacteria and Archaea)
wrp = subset_taxa(wrp,
                  domain!="Unassigned"&
                    order!="Chloroplast" &
                    order!="Mitochondria" &
                    family!="Mitochondria" &
                    domain!="Eukaryota")

# 2. Removing samples with low number of reads
# Adding the sample sums of the reads to meta data column
wrp@sam_data$sample_sums_unfiltered = as.numeric(sample_sums(wrp))
### We see a relatively continuous increase in sample, so will use 800 as the cutoff.
# Removing all samples with less than 800 reads
wrp.high <- prune_samples(sample_sums(wrp) >= 800, wrp)
# Decided to use 5000, big jump from lowest value to a pretty continuous series of values

## Making a new object with only the samples  lost (less than 5000 reads)
wrp.below800 <- prune_samples(sample_sums(wrp) < 800, wrp)

## Getting the metadata out of phyloseq for low reads obj
wrp.below800 = as.matrix(wrp.below800@sam_data)

## write file to know which samples were lost here. This is important for the methods section.
write.csv(wrp.below800, "data/processed/wrp_samples_less_than_800.csv")

## removing ASVs that are low frequency ####
### extracting otu dataframe (asv table) from phyloseq object
otutab <- as.data.frame(t(as.matrix(otu_table(wrp.high@otu_table))))

### Calcuting the sum of each row in the otuta
otutab$asv_abundance = rowSums(otutab)

### Finding the minimum value of asv_abundance 
# min(otutab$asv_abundance) # prints minimum value of ASV abundance

### Remvoing ASVs with low Frequency
#### Using 100 as threshold
otu.pruned = subset(otutab, otutab$asv_abundance>=100)

### Confirming the new minimum ASV abundance value is at the threshold
# min(otu.pruned$asv_abundance) # prints minimum value of ASV abundance

#### Removing ASV_abundance column
widthotu = ncol(otu.pruned) # finding width
## keep everything in the otu.pruned dataset except the last columns
otu.pruned = otu.pruned[,-c(widthotu)]


## removing infrequent ASVs over samples ####
### Creating a function that counts the number of occurence along rows where the number in the cell (ASV occurence) 
ASVoccur = function(x){return(sum(x>0))}

### Calculating the occurrence of each ASV in your dataframe
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)
### Investigating what the ASV occurance counts look like
# summary(otu.pruned$asv_occur_count) # prints summary

### Removing ASVs found two or less times
otu.highfreq = subset(otu.pruned, otu.pruned$asv_occur_count>2)

### Confirming filtering worked
summary(otu.highfreq$asv_occur_count)

### Removing the asv_occur_count column
otu.highfreq = otu.highfreq[,-c(widthotu)]

## data de-noising ####
otu.clean <- mutate_all(otu.highfreq, funs(ifelse(. < 3, 0, .)))

## Making a new phyloseq object of cleaned king data
cleanwrp = phyloseq(sample_data(wrp.high@sam_data),
                    tax_table(wrp.high@tax_table),
                    otu_table(as.matrix(otu.clean), taxa_are_rows = TRUE))


# creating the dephyloseq function ####
### breaks down the phyloseq object into more manageable data
dephyloseq = function(phylo_obj){
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  ## how many metadata columns you have
  metacols = ncol(meta)+1
  ## get out the otu table
  ## if your metadta is empty after running this, you need to use
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  ## get out the taxonomy file
  tax = as.data.frame(phylo_obj@tax_table)
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_name")
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_name", values_to="asv_abundance")
  ## Join the metadata and otu table with the taoxnomy table
  mot = full_join(mo, tax)
  ## Specify the output for the dephyloseq function
  output = mot
}

# choosing to use Rarefied data going forward
# leaving this here as a method to quickly decide to instead use rarefied if that is preferable
#working.wrp = rare.wrp 
working.wrp = cleanwrp

## summarize at rank 6 (or change this to be the kans you have in your dataset)
working.wrp = tax_glom(working.wrp, taxrank = "genus")
## calculate the number of reads in each sample. This is important for relative abundance calculations later
working.wrp@sam_data$read_depth = sample_sums(working.wrp)

## get working.wrp data out of phyloseq and into a dataframe
wrp.processed.df = dephyloseq(working.wrp)

# Finding lowest taxanomic data ####
# Creating function to remove blank strings and replace with NAs
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}

# Applying the above function to our dataset
wrp.processed.df <- wrp.processed.df %>%
  mutate(across(everything(), replace_empty_with_na))

# Creating Final Dataset ####
community_data <- wrp.processed.df %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(genus), genus,
                               if_else(!is.na(family), family,
                                       if_else(!is.na(order), order,
                                               if_else(!is.na(class), class,
                                                       if_else(!is.na(phylum), phylum, NA_character_))))))

# Writing to a file ####
write.csv(community_data, "data/processed/community_dataframe.csv")



### Setting up the functional trait dataset

## Importing Data ####
traits <- read.csv("data/raw/msystems.01422-21-s0005.csv")

## Condensing categories of traits into single columns ####

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

## Condensing categories into one single nitrogen-related column ####

traits <- traits %>%
  mutate_at(vars(c("dissimilatory_nitrate_reduction","assimilatory_nitrate_reduction","denitrification","nitrogen_fixation","nitrification")), ~coalesce(., "NA")) %>%
  mutate(nitrogen_cycling = if_else(rowSums(select(., c("dissimilatory_nitrate_reduction","assimilatory_nitrate_reduction","denitrification","nitrogen_fixation","nitrification")) == "Y", na.rm = TRUE) > 0, "Y", ""))

## Creating a new condensed dataframe, with only one variable for nitrogen cycling ####

condensed_traits <- traits[1:66,c("MAG_name", "phylum","class","order","family","genus","nitrogen_cycling")]

## Finding lowest taxanomic data ####
# Creating function to remove blank strings and replace with NAs
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}

replace_na_with_n <- function(x) {
  x[is.na(x)] <- "N"
  return(x)
}

# Applying the above function to our dataset
condensed_traits <- condensed_traits %>%
  mutate(across(everything(), replace_empty_with_na))

trait_data <- condensed_traits %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(genus), genus,
                               if_else(!is.na(family), family,
                                       if_else(!is.na(order), order,
                                               if_else(!is.na(class), class,
                                                       if_else(!is.na(phylum), phylum, NA_character_))))))

trait_data <- trait_data %>%
  mutate(nitrogen_cycling = replace_na_with_n(nitrogen_cycling))

# Dealing with repeat cases of lowest_rank
# going to merge these observations into a single value for each lowest_rank
# MAJOR DEBATE HERE: Choosing to make merged rows, such that if the rows ever say Y to N cycling, they always do

trait_data <- trait_data %>% arrange(desc(nitrogen_cycling))
trait_data <- trait_data %>% 
  group_by(lowest_rank) %>% 
  slice(1) %>% 
  ungroup()


### Merging the two datasets

# trait_data is the dataframe with known taxa and traits
# community_data is the dataframe is the collection of data relating to commmunity composition and shifts based on treatment.

# We will want to join these two dataframes using some version of the join() function
# inner_join() keeps only observations that exist in both - might lose a lot of data
# left_join() keeps observations in the left (first) dataframe
# full_join() keeps all observations in both

# My (elias') thought: either inner_join - and work with only the ones that have both
# or, we could try left join, and keep all the observations from the original data, and throw out non-relevant trait_data

all_merged <- left_join(community_data, trait_data, by = "lowest_rank")

# Creating two sets of data one for meristem v tip and one for host type
# In either case filtering will need to happen
## Filtering for relevant samples
# In species_data - only want macro and nereo (maybe control)
# In blade_data - only want nereo base and nereo tip (whats the differnec between "tip" and "tip_swab")
## Filtering for relevant sites 
# Either question can only be asked using data from specific sites - need to filter for these locations

both_species_locations <- c("Bullman", "Cape_Johnson", "Destruction_Island", 
                            "Koitlah", "Sekiu")
base_and_tip <- c("Nereocystis_tip", "Nereocystis_base", "Nereocystis_tip_swab",
                  "Nereocystis_base_swab")

# The below makes species_data, with only samples from the above listed locations, where both kelp species were sampled
species_data <- all_merged %>%
  select(Row.names, description, host_common_name, location, sample_type_user, asv_abundance, lowest_rank, nitrogen_cycling) %>%
  filter(description == "Macrocystis" | description == "Nereocystis") %>%
  filter(location %in% both_species_locations)

# The below makes blade_data, with only samples from Tatoosh and from the blade meristem and tip, also makes a new simplified 
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


