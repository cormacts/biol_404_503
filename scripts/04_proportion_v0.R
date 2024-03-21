####
# BIOL 403/503 project
# Lydia
# Working on adding proportions to the analysis

### Species comparison

# Grouping by sample, so we have multiple samples per species of kelp
new_sum_sp_data <- species_data %>%
  group_by(Row.names, nitrogen_cycling,description) %>%
  summarise_at(vars(asv_abundance),
               list(sum_abundance = sum))
# Calculating total abundance
proportion_data <- new_sum_sp_data %>%
  group_by(Row.names,description) %>%
  summarise_at(vars(sum_abundance),
               list(total_abundance = sum))

combined_proportions <- merge(new_sum_sp_data,proportion_data)

# Next step - divide sum_abundance by total_abundance for each sample
               
               