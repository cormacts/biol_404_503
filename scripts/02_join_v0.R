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
# or, we could try left join, and keep all the observations from the original data, and throw out non-relevant trat_data
