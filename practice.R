# RD Diss data analysis

# Workflows:
# Clean and aggregate aphid data by plot - DONE
# Clean and aggregate parasitoid wasp data - DONE
# Merge dataasets - DONE
# Prepare wildflower strip data:
      # Want to end up with a single metirc for diversity as a proxy value for nectar sources - IN PROGRESS
      # Use distance from wildflower strip as an effect in model - TODO
# Aggregate in-plot wildflower data:
      # Assess homogeneity of wildflowers across the plots - IN PROGRESS
      # Use Wildflower variability across plots as a predictor in the model - precise metric to be determined
# Example model: glm(mean_infested ~ mean_colemani + wsf_metric + inplot_wf, data = merged_data, family = "nb" )

#rm(list=ls())
# load packages
library(tidyverse)
library(glmmTMB) # models
library(lme4) # models?
library(MuMIn) # AICc
library(readxl)
library(vegan)  # For diversity indices?
 
#read in csv's
wasp_data <- read_excel("data/RMD_sticky_trap_data.xlsx")
aphid_data <- read_excel("data/aphid_data.xlsx")
wfs_data <- read_excel("data/Berryhill_Wildflower_Survey_copy.xlsx", sheet = "wfs_data")
plot_wf_data <- read_excel("data/Berryhill_Wildflower_Survey_copy.xlsx", sheet = "Wheeling")

# explore data
view(wasp_data)
view(aphid_data)
view(wfs_data)
view(plot_wf_data)
str(wasp_data)
str(aphid_data)
str(plot_wf_data)
# cleaning ----
wasp_data <- wasp_data %>% 
  mutate(across(c(set_no, n_treatment, crop, margin, con_rel), as.factor))
aphid_data <- aphid_data %>%
  mutate(across(c(set_no,sample_point), as.factor))

mean(wasp_data$aphidius_colemani)   #  mean
var(wasp_data$aphidius_colemani) # variance (negative binomial check)

# average aphid and wasp abundance across time (aggregate by plot) ----

aphid_avg <- aphid_data %>% 
      group_by(plot_id) %>% 
  summarise(mean_aphid = mean(no_infested))
#view(aphid_avg)
wasp_avg <- wasp_data %>% 
      group_by(plot_id) %>%
  summarise(
    mean_ervi = mean(aphidius_ervi),
    mean_colemani = mean(aphidius_colemani),
    mean_ichneumonid = mean(ichneumonid_spp.),
    mean_other = mean(other_parasitica),
    .groups = "drop"
  )
#view(wasp_avg)
aphid_aggregated <- aphid_data %>%
  group_by(plot_id) %>%
  summarise(
    total_infested = sum(no_infested),  # Total across all time/sample points
    mean_infested = mean(no_infested),  # Mean per sample point
    .groups = "drop"
  )

wasp_aggregated <- wasp_data %>%
  group_by(plot_id) %>%
  summarise(                                 # New data frame with average abundance of each species
    mean_ervi = sum(aphidius_ervi),
    mean_colemani = sum(aphidius_colemani),
    mean_ichneumonid = sum(ichneumonid_spp.),
    mean_other = sum(other_parasitica),
    .groups = "drop"                         # not quite sure why I use drop here lol
  ) %>%
  # Convert to numeric and remove "P" from plot_id data
  mutate(plot_id = as.numeric(gsub("P", "", plot_id))) %>%
  arrange(plot_id)

# Combine data sets
merged_data <- inner_join(
  aphid_aggregated, 
  wasp_aggregated, 
  by = "plot_id"
)
view(wasp_aggregated)
#str(wasp_aggregated)
view(merged_data)
# wildflower strip data wrangling ----

# Convert Q1-Q10 columns to character
wfs_data <- wfs_data %>%
  mutate(across(starts_with("Q") & !ends_with("m_c_r"), as.character)) %>%
  mutate(across())

# create new column of quadrat braun-Blanquet scores
wildflower_bb <- wfs_data %>%
   select(`common name`, Species, starts_with("Q") & !ends_with("m_c_r")) %>%
   pivot_longer(cols = -c(`common name`, Species),
   names_to = "quadrat",
   values_to = "braun_blanquet")

# so the same for quadrat mid-point cover values

wildflower_mc <- wfs_data %>%
  select(`common name`, Species, ends_with("m_c_r")) %>%
  pivot_longer(cols = -c(`common name`, Species),
    names_to = "quadrat",
    values_to = "midpoint_cover"
    ) %>%
  mutate(quadrat = gsub("_m_c_r", "", quadrat))


# merge the two datasets
wfs_long <- wildflower_bb %>%
    left_join(wildflower_mc, by = c("common name", "Species", "quadrat")) %>%
  arrange(quadrat)
view(wfs_long)
#str(wfs_long)

# collate floral cover per wfs quadrat

quadrat_metrics <- wfs_long %>%
  group_by(quadrat) %>%
  summarise(
    total_cover = sum(midpoint_cover, na.rm = TRUE),  # Total cover per quadrat
    shannon_div = diversity(midpoint_cover, index = "shannon")  # Diversity per quadrat
  )
view(quadrat_metrics)
# in-plot wildflower data wrangling ----

plot_wf_data <- plot_wf_data %>%
  mutate(across(starts_with("P"), as.character)) %>%
  pivot_longer(cols = starts_with("P"),
    names_to = "plot_position",
    values_to = "braun_blanquet") %>%
  separate(col = plot_position,
    into = c("plot_id", "quadrat"),  # Quoted column names
    sep = "C") %>%                        # Get rid of text from  "C" (e.g., "P1C2" → "P1")
  mutate(
    plot_id = as.numeric(gsub("P", "", plot_id)),  # Convert "P1" → 1
    quadrat = paste0("C", quadrat)                # Restore quadrat label (e.g., "C2")
  ) %>%

  mutate(
    midpoint_cover = case_when(
      braun_blanquet == "+" ~ 0.1,
      braun_blanquet == "1" ~ 2.5,
      braun_blanquet == "2" ~ 15,
      braun_blanquet == "3" ~ 37.5,
      braun_blanquet == "4" ~ 62.5,
      braun_blanquet == "5" ~ 87.5,
      TRUE ~ 0
    ))# %>%
  mutate(mc_scaled = midpoint_cover * (1 / 0.15))
  
view(plot_wf_data)
str(plot_wf_data)

plot_metrics <- plot_wf_data %>%
  group_by(plot_id) %>%
  summarise(
    total_cover = sum(mc_scaled, na.rm = TRUE),
    shannon_div = vegan::diversity(mc_scaled, index = "shannon")
  )


view(plot_metrics)
# Tidy, scale, and aggregate
plot_metrics <- plot_wf_data %>%
  mutate(across(starts_with("P"), as.character)) %>%
  pivot_longer(cols = starts_with("P"), names_to = "plot_position", values_to = "braun_blanquet") %>%
  separate(col = plot_position, into = c("plot_id", "quadrat"), sep = "C") %>%
  mutate(
    plot_id = as.numeric(gsub("P", "", plot_id)),
    quadrat = paste0("C", quadrat),
    midpoint_cover = case_when(
      braun_blanquet == "+" ~ 0.1,
      braun_blanquet == "1" ~ 2.5,
      braun_blanquet == "2" ~ 15,
      braun_blanquet == "3" ~ 37.5,
      braun_blanquet == "4" ~ 62.5,
      braun_blanquet == "5" ~ 87.5,
      TRUE ~ 0
    ),
    midpoint_cover_scaled = midpoint_cover * (1 / 0.15)
  ) %>%
  group_by(plot_id) %>%
  summarise(
    total_cover = sum(midpoint_cover_scaled, na.rm = TRUE),
    shannon_div = vegan::diversity(midpoint_cover, index = "shannon")
  )
