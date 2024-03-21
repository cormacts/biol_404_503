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

# Adding columns calculating 
# 1. Total abundance of all ASVs that are taxa covered by Weigel 2022 (functional traits paper)
# 2. Total abundance of all ASVs covered by Weigel 2019 (original paper)
# 3. Proportions for Y, N, and Unknown nitrogen cycling traits, out of either the 2022 or 2019 total
pivoted_sum_sp_data$total_abundance_Weigel2022 = pivoted_sum_sp_data$Yes+pivoted_sum_sp_data$No
pivoted_sum_sp_data$total_abundance_Weigel2019 = pivoted_sum_sp_data$Yes+pivoted_sum_sp_data$No+pivoted_sum_sp_data$Unknown
pivoted_sum_sp_data$proportion_Y_22 = pivoted_sum_sp_data$Yes/pivoted_sum_sp_data$total_abundance_Weigel2022
pivoted_sum_sp_data$proportion_N_22 = pivoted_sum_sp_data$No/pivoted_sum_sp_data$total_abundance_Weigel2022
pivoted_sum_sp_data$proportion_Y_19 = pivoted_sum_sp_data$Yes/pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_sp_data$proportion_N_19 = pivoted_sum_sp_data$No/pivoted_sum_sp_data$total_abundance_Weigel2019
pivoted_sum_sp_data$proportion_NA_19 = pivoted_sum_sp_data$Unknown/pivoted_sum_sp_data$total_abundance_Weigel2019

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
# 2. Total abundance of all ASVs covered by Weigel 2019 (original paper)
# 3. Proportions for Y, N, and Unknown nitrogen cycling traits, out of either the 2022 or 2019 total
pivoted_sum_blade_data$total_abundance_Weigel2022 = pivoted_sum_blade_data$Yes+pivoted_sum_blade_data$No
pivoted_sum_blade_data$total_abundance_Weigel2019 = pivoted_sum_blade_data$Yes+pivoted_sum_blade_data$No+pivoted_sum_blade_data$Unknown
pivoted_sum_blade_data$proportion_Y_22 = pivoted_sum_blade_data$Yes/pivoted_sum_blade_data$total_abundance_Weigel2022
pivoted_sum_blade_data$proportion_N_22 = pivoted_sum_blade_data$No/pivoted_sum_blade_data$total_abundance_Weigel2022
pivoted_sum_blade_data$proportion_Y_19 = pivoted_sum_blade_data$Yes/pivoted_sum_blade_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_N_19 = pivoted_sum_blade_data$No/pivoted_sum_blade_data$total_abundance_Weigel2019
pivoted_sum_blade_data$proportion_NA_19 = pivoted_sum_blade_data$Unknown/pivoted_sum_blade_data$total_abundance_Weigel2019

## Making a new sample ID column with species in it
pivoted_sum_sp_data$sample_ID = paste(pivoted_sum_sp_data$description,pivoted_sum_sp_data$Row.names)

## Making dataframe to use for plotting
plot_sp_proportions <- pivoted_sum_sp_data[,c("sample_ID","proportion_Y_22","proportion_N_22")]
plot_sp_proportions <- plot_sp_proportions %>%
  pivot_longer(!sample_ID, names_to = "nitrogen_cycling",values_to = "proportions_22")
plot_sp_proportions <- plot_sp_proportions %>% 
  separate(sample_ID,c("species", "sample"), sep = " ", remove = TRUE)

## Plotting code from Annie + taxaplots code from Lab 5:
ggplot(plot_sp_proportions, aes(x=sample, y=proportions_22,
                  fill=nitrogen_cycling))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("lightblue", "yellow1"))+
  guides(fill=guide_legend(ncol=2))+
  facet_grid(.~species, scales="free", space="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size=10, face="bold"),
        strip.text = element_text(color="black", size=10),
        legend.text=element_text(size=6),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank())+
  labs(y="Relative abundance", x="Sample", fill="Nitrogen cycling")
