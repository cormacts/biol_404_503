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

