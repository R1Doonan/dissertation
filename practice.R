# RMD
# Diss data analysis
 rm(list=ls())
# load packages
library(tidyverse)
library(glmmTMB) # models
library(MuMIn) # AICc
library(readxl)
library(vegan)  # For diversity indices
 
#read in csv's
wasp_data <- read_excel("data/RMD_sticky_trap_data.xlsx")
aphid_data <- read_excel("data/aphid_data.xlsx")
wfs_data <- read_excel("data/Berryhill_Wildflower_Survey_copy.xlsx", sheet = "wfs_data")
plot_wf_data <- read_excel("data/Berryhill_Wildflower_Survey_copy.xlsx", sheet = "Wheeling-plotcombined")

# explore data
view(wasp_data)
view(aphid_data)
view(wfs_data)
view(plot_wf_data)
str(wasp_data)
str(aphid_data)
str(plot_wf_data)
# cleaning (mutate to covert to factors)
wasp_data <- wasp_data %>% 
  mutate(across(c(set_no, n_treatment, crop, margin, con_rel), as.factor))
aphid_data <- aphid_data %>%
  mutate(across(c(set_no,sample_point), as.factor))

mean(wasp_data$aphidius_colemani)   #  mean
var(wasp_data$aphidius_colemani) # variance (negative binomial check)

# average aphid and wasp abundance across time ----

aphid_avg <- aphid_data %>% 
      group_by(plot_id) %>% 
  summarise(mean_aphid = mean(no_infested))
view(aphid_avg)
wasp_avg <- wasp_data %>% 
      group_by(plot_id) %>%
  summarise(
    mean_ervi = mean(aphidius_ervi),
    mean_colemani = mean(aphidius_colemani),
    mean_ichneumonid = mean(ichneumonid_spp.),
    mean_other = mean(other_parasitica),
    .groups = "drop"
  )
view(wasp_avg)
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
  arrange(plot_id)
  #select(-plot_id)  # Remove original ID column
merged_data <- inner_join(
  aphid_aggregated, 
  wasp_aggregated, 
  by = "plot_id"
)
view(wasp_aggregated)
str(wasp_aggregated)
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
    left_join(wildflower_mc, by = c("common name", "Species", "quadrat"))
view(wfs_long)
str(wfs_long)

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
  mutate(across(starts_with("p"), as.character)) %>%
  mutate(across()) %>%
  pivot_longer(cols = starts_with("P"),
    names_to = "plot_position",
    values_to = "braun_blanquet")

view(plot_wf_data)


