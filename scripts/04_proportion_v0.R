####
# BIOL 403/503 project
# Lydia
# Working on adding proportions to the analysis

### Species comparison

# Grouping by sample (assuming that it is being stored in Row.names) - so multiple samples per species of kelp
new_sum_sp_data <- species_data %>%
  group_by(Row.names, nitrogen_cycling,description) %>%
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

### Making dataframe to use for plotting
pivoted_sum_sp_data$sample_ID = paste(pivoted_sum_sp_data$description,pivoted_sum_sp_data$Row.names)
plot_sp_proportions <- pivoted_sum_sp_data[,c("sample_ID","proportion_Y_22","proportion_N_22")]
plot_sp_proportions <- plot_sp_proportions %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling",values_to = "proportions_22")
plot_sp_proportions <- plot_sp_proportions %>% 
  separate(sample_ID,c("species", "sample"), sep = " ", remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
#png(file = "figures/sp_proportion_taxaplot.png")
png(
  file = "figures/sp_proportion_taxaplot.png",
  width     = 5,
  height    = 5,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
ggplot(plot_sp_proportions, aes(x=sample, y=proportions_22,
                                fill=nitrogen_cycling))+
  geom_bar(stat = "identity")+
  scale_fill_manual(labels = c("No", "Yes"), values = c("lightblue", "yellow1"))+
  guides(fill=guide_legend(ncol=2))+
  facet_grid(.~species, scales="free", space="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=10),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank())+
  labs(y="Relative abundance", x="Sample", fill="Nitrogen cycling")
dev.off()


### Blade location / successional stage comparison

# Grouping by sample (assuming that it is being stored in Row.names) - so multiple samples per blade location
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

## Making dataframe to use for plotting
pivoted_sum_blade_data$sample_ID = paste(pivoted_sum_blade_data$sample_type,pivoted_sum_blade_data$Row.names)
plot_blade_proportions <- pivoted_sum_blade_data[,c("sample_ID","proportion_Y_22","proportion_N_22")]
plot_blade_proportions <- plot_blade_proportions %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling",values_to = "proportions_22")
plot_blade_proportions <- plot_blade_proportions %>% 
  separate(sample_ID,c("blade_location", "sample"), sep = " ", remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
#png(file = "figures/b_proportion_taxaplot.png")
png(
  file = "figures/b_proportion_taxaplot.png",
  width     = 5,
  height    = 5,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
ggplot(plot_blade_proportions, aes(x=sample, y=proportions_22,
                                   fill=nitrogen_cycling))+
  geom_bar(stat = "identity")+
  scale_fill_manual(labels = c("No", "Yes"), values = c("lightblue", "yellow1"))+
  guides(fill=guide_legend(ncol=2))+
  facet_grid(.~blade_location, scales="free", space="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=10),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank())+
  labs(y="Relative abundance", x="Sample", fill="Nitrogen cycling")
dev.off()

##### Same procedure but with full Weigel 2019 dataset and taxa with unknown nitrogen-cycling traits

### Species comparison

## Making dataframe to use for plotting
plot_sp_proportions_2 <- pivoted_sum_sp_data[,c("sample_ID","proportion_Y_19","proportion_N_19","proportion_NA_19")]
plot_sp_proportions_2 <- plot_sp_proportions_2 %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling",values_to = "proportions_19")
plot_sp_proportions_2 <- plot_sp_proportions_2 %>% 
  separate(sample_ID,c("species", "sample"), sep = " ", remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
##png(file = "figures/proportion_NAs_taxaplot.png")

png(
  file = "figures/proportion_NAS_taxaplot.png",
  width     = 5,
  height    = 5,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)

ggplot(plot_sp_proportions_2, aes(x=sample, y=proportions_19,
                                fill=nitrogen_cycling))+
  geom_bar(stat = "identity")+
  scale_fill_manual(labels = c("No", "Unknown", "Yes"), values = c("lightblue","violet","yellow1"))+
  guides(fill=guide_legend(ncol=2))+
  facet_grid(.~species, scales="free", space="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size=10, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=10),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank())+
  labs(y="Relative abundance", x="Sample", fill="Nitrogen cycling")
dev.off()

#### Blade location comparison

## Making dataframe to use for plotting
plot_blade_proportions_2 <- pivoted_sum_blade_data[,c("sample_ID","proportion_Y_19","proportion_N_19","proportion_NA_19")]
plot_blade_proportions_2 <- plot_blade_proportions_2 %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling",values_to = "proportions_19")
plot_blade_proportions_2 <- plot_blade_proportions_2 %>% 
  separate(sample_ID,c("blade_location", "sample"), sep = " ", remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
#png(file = "figures/b_proportion_NAs_taxaplot.png")

png(
  file = "figures/b_proportion_NAS_taxaplot.png",
  width     = 5,
  height    = 5,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
ggplot(plot_blade_proportions_2, aes(x=sample, y=proportions_19,
                                   fill=nitrogen_cycling))+
  geom_bar(stat = "identity")+
  scale_fill_manual(labels = c("No", "Unknown", "Yes"), values = c("lightblue","violet","yellow1"))+
  guides(fill=guide_legend(ncol=2))+
  facet_grid(.~blade_location, scales="free", space="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size=10, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=10),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank())+
  labs(y="Relative abundance", x="Sample", fill="Nitrogen cycling")
dev.off()
