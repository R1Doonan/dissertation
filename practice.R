# RMD
# Diss data analysis
 rm(list=ls())
# load packages
library(tidyverse)
library(glmmTMB)
library(MuMIn)
#read in csv's
wasp_data <- read_excel("data/RMD_sticky_trap_data.xlsx")
aphid_data <- read_excel("data/aphid_data.xlsx")

# explore data
view(wasp_data)
view(aphid_data)
str(wasp_data)
str(aphid_data)

# cleaning (mutate)
wasp_data <- wasp_data %>% 
  mutate(across(c(set_no, n_treatment, crop, margin, con_rel), as.factor))
aphid_data <- aphid_data %>%
  mutate(across(c(set_no,sample_point), as.factor))

mean(wasp_data$aphidius_colemani)   #  mean
var(wasp_data$aphidius_colemani) # variance



# average aphid and wasp data across time ----

aphid_avg <- aphid_data %>% 
      group_by(plot_id) %>% 
  summarise(mean_aphid = mean(no_infested))

wasp_avg <- wasp_data %>% 
      group_by(plot_id) %>%
  summarise(
    mean_ervi = mean(aphidius_ervi),
    mean_colemani = mean(aphidius_colemani),
    mean_ichneumonid = mean(ichneumonid_spp.),
    mean_other = mean(other_parasitica),
    .groups = "drop"
  )

aphid_aggregated <- aphid_data %>%
  group_by(plot_id) %>%
  summarise(
    total_infested = sum(no_infested),  # Total across all time/sample points
    mean_infested = mean(no_infested),  # Mean per sample point
    .groups = "drop"
  )

wasp_aggregated <- wasp_data %>%
  group_by(plot_id) %>%
  summarise(
    mean_ervi = mean(aphidius_ervi),
    mean_colemani = mean(aphidius_colemani),
    mean_ichneumonid = mean(ichneumonid_spp.),
    mean_other = mean(other_parasitica),
    .groups = "drop"
  ) %>%
  # Convert plot_id "P1" to numeric plot_id 1
  mutate(plot_id = as.numeric(gsub("P", "", plot_id))) %>%
  arrange(plot_id) %>%
  #select(-plot_id)  # Remove original ID column
merged_data <- inner_join(
  aphid_aggregated, 
  wasp_aggregated, 
  by = "plot_id"
)
view(wasp_aggregated)
str(wasp_aggregated)

