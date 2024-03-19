#############
# BIOL 403 503 Project
# Script: 00_setup_v0
# By: Elias Bowman
# 2024-02-12
# 1st script of BIOL 403 Lab project, working with project data set
# Begin by importing relevant packages 
# Then filter the data set: 
#   Remove off target data, failed samples, wroing ASVs, few samples, and low ASV counts
###########

# Initial Setup ####
## Loading Libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(plyr)
library(qualpalr)
library(ggpubr)

## Reading in data from mar2018
wrp = readRDS("data/raw/W2019_and_RP2022_unfiltered_phyloseq.RDS")

# Filtering ####
## remove off target taxa ####
### We want only Bacteria and Archea within our data
wrp = subset_taxa(wrp,
                  domain!="Unassigned"&
                    order!="Chloroplast" &
                    order!="Mitochondria" &
                    family!="Mitochondria" &
                    domain!="Eukaryota")

## removing samples with low number of reads ####
# sampleSums(wrp) # This code prints the sample sums of 

### Add the sample sums of the reads to meta data column
wrp@sam_data$sample_sums_unfiltered = as.numeric(sample_sums(wrp))
#### relatively continuous increase in sample, going to use 800 as the cutoff 

## Remove all samples with less than 800 reads
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
otu.pruned = otu.pruned[,-c(widthotu)]


## removing infrequent ASVs over samples ####
### Creating a function that counts the number of occurence along rows where the number in the cell (ASV occurence) 
ASVoccur = function(x){return(sum(x>0))}

### Calculating the occurrence of each ASV in your dataframe
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)
### Investigating what the ASV occurance counts look like
summary(otu.pruned$asv_occur_count) # prints summary

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

## phyloseq-class experiment-level object
## otu_table() OTU Table: [ 337 taxa and 248 samples ]
## sample_data() Sample Data: [ 248 samples by 31 sample variables ]
## tax_table() Taxonomy Table: [ 337 taxa by 8 taxonomic ranks ]

## save rarefied dataframe
write_rds(cleanwrp, "data/processed/wrp_filtered.RDS")

# Creating Taxaplots ####
## creating the dephyloseq function
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

## summarize at rank 6 (or change this to be the kans you have in your dataset)
cleanwrp = tax_glom(cleanwrp, taxrank = "genus")
## calculate the number of reads in each sample. This is important for relative abundance calculations later
cleanwrp@sam_data$read_depth = sample_sums(cleanwrp)

## get cleanwrp data out of phyloseq and into a dataframe
cleanwrp.df = dephyloseq(cleanwrp)
write.csv(cleanwrp.df, "data/processed/cleanwrp_dataframe.csv")


## Finding lowest taxanomic data ####
# Creating function to remove blank strings and replace with NAs
replace_empty_with_na <- function(x) {
  x[x == ""] <- NA
  return(x)
}

# Applying the above function to our dataset
cleanwrp.df <- cleanwrp.df %>%
  mutate(across(everything(), replace_empty_with_na))

community_data <- cleanwrp.df %>%
  rowwise() %>%
  mutate(lowest_rank = if_else(!is.na(genus), genus,
                       if_else(!is.na(family), family,
                       if_else(!is.na(order), order,
                       if_else(!is.na(class), class,
                       if_else(!is.na(phylum), phylum, NA_character_))))))
