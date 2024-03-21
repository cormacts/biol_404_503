#############
# BIOL 403 503 Project
# Script: 02_join_v0
# By: Lydia, Elias
# 2024-03-18
# This script will join the functional trait data set and the community data set
#############

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


