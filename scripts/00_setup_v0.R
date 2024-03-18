#############
# BIOL 403 503 Project
# Script: 00_setup_v0
# By: Elias Bowman
# 2024-02-12
# 1st script of BIOL 403 Lab project, working with project data set
# Begin by importing relevant packages 
# Then filter the data set: 
#   Remove off target data, failed samples, wroing ASVs, few samples, and low ASV counts
#############

# Initial Setup ####
## Loading Libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(plyr)
library(qualpalr)
library(ggpubr)

## Setting my Working Directory
# setwd("C:/Users/elias/OneDrive/Documents/University/BIOL 403/BIOL403_proj")

## Reading in data from mar2018
wrp = readRDS("data/W2019_and_RP2022_unfiltered_phyloseq.RDS")

## Viewing the 3 parts of the phyloseq object
# metadata table
# View(wrp@sam_data)
# view taxonomy table
# View(wrp@tax_table)
# otu table
# View(wrp@otu_table)

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
sampleSums(wrp)
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
write.csv(wrp.below800, "wrp_samples_less_than_800.csv")

## removing ASVs that are low frequency ####
### extracting otu dataframe (asv table) from phyloseq object
otutab <- as.data.frame(t(as.matrix(otu_table(wrp.high@otu_table))))

### Calcuting the sum of each row in the otuta
otutab$asv_abundance = rowSums(otutab)

### Finding the minimum value of asv_abundance 
min(otutab$asv_abundance)

### Remvoing ASVs with low Frequency
#### Using 100 as threshold
otu.pruned = subset(otutab, otutab$asv_abundance>=100)

### Confirming the new minimum ASV abundance value is at the threshold
min(otu.pruned$asv_abundance)

#### Removing ASV_abundance column
widthotu = ncol(otu.pruned) # finding width
otu.pruned = otu.pruned[,-c(widthotu)]


## removing infrequent ASVs over samples ####
### Creating a function that counts the number of occurence along rows where the number in the cell (ASV occurence) 
ASVoccur = function(x){return(sum(x>0))}

### Calculating the occurrence of each ASV in your dataframe
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)
### Investigating what the ASV occurance counts look like
summary(otu.pruned$asv_occur_count)

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
write_rds(cleanwrp, "wrp_filtered.RDS")

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
  ## Specify the output for the dephyloseq funciton
  output = mot
}

## summarize at rank 6 (or change this to be the kans you have in your dataset)
cleanwrp = tax_glom(cleanwrp, taxrank = "genus")
## calculate the number of reads in each sample. This is important for relative abundance calculations later
cleanwrp@sam_data$read_depth = sample_sums(cleanwrp)

## get cleanwrp data out of phyloseq and into a dataframe
cleanwrp.df = dephyloseq(cleanwrp)


# Taxaplots ####
# normalize the # of reads in each sample by rel. abundance to count ASV
## caluclate relative abundance of each rank 6/genus within each sample
cleanwrp.df$relativeabundance = as.numeric(cleanwrp.df$asv_abundance)/as.numeric(cleanwrp.df$read_depth)

# making names for the plots
cleanwrp.df$plotnames = paste0(cleanwrp.df$order, ";", cleanwrp.df$genus)



# creating loop to make plots ####
## summarize data by substrate type. This will let you make one taxaplot for each substrate, so 3 plots ## you will need to change substrate to YOUR variable of interest
wrp.sum = ddply(cleanwrp.df, c("location", "plotnames"),
                summarise,
                sum = sum(relativeabundance))

## get list of the substrate types. This is what the loop will use
samplegroups = unique(wrp.sum$location)

## sort data by relative abundance. This is how the loop will pick the most abundant taxa
sorted = wrp.sum[order(-wrp.sum$sum),]

## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in samplegroups) {
  for(j in i) {
    ## subset dataframe by samples
    #!# Remeber to change the substrate to your group!
    sample = subset(sorted, sorted$location %in% c(j))
    ## get top 15 genera
    top = sample[c(1:15),]
    ## save list of top abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    ## close loop
  }
}

## add identifier for top 15 taxa
top.df$place = "top_15"

# Combine your top 15 taxa for each substrate type (from the loop) to your entire dataset
## join the top taxa and existing dataframe
alldata = full_join(cleanwrp.df, top.df)

## make the empty "place" cells say bottom. This workes because we used full_join
alldata$place = replace(alldata$place, is.na(alldata$place), "bottom")
## replace plot_names that have bottom taxa as their "place" with Other
alldata[alldata$place == "bottom",]$plotnames <- "Others"

# Picking out colours ####
# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))
# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(15)
# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")
# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")
# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")
# 6. Join color list and taxa names
topcolors = cbind(tops, hex)
# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"

# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

# Order Taxaplot ####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$relativeabundance),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))


# Making the Taxaplot ####
## these plots get pretty big, so let's only plot the macrocystis samples
macro = subset(alldata, alldata$description =="Macrocystis")
## make the plot
ggplot(macro, aes(x=as.character(Row.names), y=as.numeric(relativeabundance),
                  fill=as.factor(plotnames)))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=plotcolors)+
  guides(fill=guide_legend(ncol=2))+
  facet_grid(.~location, scales="free", space="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size=10, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=6),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank())+
  labs(y="Relative Abundance", x="Sample", fill="Taxa")
