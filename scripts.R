# RD Diss data analysis

# Workflows:
# Clean and aggregate aphid data by plot - DONE
# Clean and aggregate parasitoid wasp data - DONE
# Merge data sets - DONE
# Prepare wildflower strip data:
      # Want to end up with a single metric for diversity as a proxy value for nectar sources - IN PROGRESS
      # Use distance from wildflower strip as an effect in model - TODO
# Aggregate in-plot wildflower data:
      # Assess homogeneity of wildflowers across the plots - IN PROGRESS 
      # Use Wildflower variability across plots as a predictor in the model - precise metric to be determined
# Example model: glm(mean_infested ~ mean_colemani + wsf_metric + inplot_wf, data = merged_data, family = "nb" )

rm(list=ls())
# load packages
library(tidyverse)
library(glmmTMB) # models
library(lme4) # models?
library(MuMIn)
library(ggeffects) # model visualisations
library(DHARMa) # model diagnostics 
library(readxl) # interpret spreadsheets
library(vegan) # For diversity indices? - seems heavily biased for this type of research

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
  mutate(plot_id = as.numeric(gsub("P", "", plot_id))) %>%   # Convert to numeric and remove "P" from plot_id data
  mutate(across(c(plot_id,set_no, time_point,n_treatment, crop, margin_distance, augmentation), as.factor))
aphid_data <- aphid_data %>%
  mutate(across(c(plot_id,set_no,time_point,sample_point), as.factor))
# augmentation effect on parasitoid abundance
str(wasp_data)
# Wrangling ----


#view(aphid_avg)
wasp_aggregated <- wasp_data %>% 
      group_by(plot_id, time_point) %>%
  arrange(plot_id)
view(wasp_aggregated)

aphid_aggregated <- aphid_data %>%
  group_by(plot_id, time_point) %>%
  summarise(
    #aphid_total_infested = sum(no_infested),  # Total across all time/sample points
    sum_infested = sum(no_infested),  # Mean per sample point
    .groups = "drop"
  )
view(aphid_aggregated)

# Combine data sets
merged_data <- inner_join(
  aphid_aggregated, 
  wasp_aggregated, 
  by = c("plot_id", "time_point")
)
#str(wasp_aggregated)
view(merged_data)
str(merged_data)
# explore data ----
check_overdispersion <- function(response_var) {
  var <- merged_data[[response_var]]
  dispersion <- var(var) / mean(var)
  cat(response_var, "dispersion ratio:", dispersion, "\n")
}

check_overdispersion("aphidius_colemani")  # 1.98>1.5 so over dispersed 
check_overdispersion("sum_infested") # 7.2 is massively over dispersed

table(merged_data$aphidius_colemani)  # Replace with your response
# plot incrop aphid vs colemani
ggplot(merged_data, aes(x = sum_infested, y = aphidius_colemani)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE)
# plot colemani vs augmentation
ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot() +
  geom_smooth(method = "loess", se = FALSE)
# plot colemani abundance vs wfs distance
ggplot(merged_data, aes(x = margin_distance, y = aphidius_colemani)) +
  geom_boxplot() +
  geom_smooth(method = "loess", se = FALSE)

# plot aphid abundance against wfs distance
ggplot(merged_data, aes(x = margin_distance, y = sum_infested)) +
  geom_boxplot() +
  geom_smooth(method = "loess", se = FALSE)

# plot aphid abundance vs augmentation
ggplot(merged_data, aes(x = augmentation, y = sum_infested)) +
  geom_boxplot() +
  geom_smooth(method = "loess", se = FALSE)

# ok model time! - for each rq lets start iterating 

# RQ1: augmentation and wildflowers effect on colemani ----
r1_mod <- glmer.nb(
  aphidius_colemani ~ augmentation * margin_distance + sum_infested + (1|time_point) +
    (1|plot_id),
  data = merged_data
)
r1_mod <- glmer.nb(
  aphidius_colemani ~ augmentation * margin_distance + sum_infested + (1|time_point),
  data = merged_data, family = "negativebinomial"
)
r1_mod <- glmer.nb(
  aphidius_colemani ~ augmentation * margin_distance + sum_infested + (1|plot_id),
  data = merged_data
)
r1_mod_final <- glmmTMB(
  aphidius_colemani ~ augmentation + margin_distance + sum_infested + (1 | plot_id),
  family = nbinom2,
  data = merged_data
)

simulateResiduals(r1_mod_final) %>% plot()
summary(r1_mod_tmb)
plot(r1_mod_final)
# RQ1 visualisations ----


# Get predicted marginal means
aug_margin_effects <- ggpredict(
  r1_mod_final,
  terms = c("augmentation", "margin_distance", "sum_infested")
)
aug_margin_effects < - aug_margin_effects %>%
  view(aug_margin_effects)
# Plot
ggplot(aug_margin_effects, aes(x = x, y = predicted, color = group)) +
  geom_boxplot(size = 0.5) +
  labs(
    x = "Augmentation",
    y = "Predicted Aphidius colemani Abundance"
  ) +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) +
  theme_minimal()
# RQ2: Augmentation and wildflower strip effect on aphid population ----
r2_mod <- glmer.nb(
  sum_infested ~ augmentation * margin_distance + sum_infested + time_point +
    (1 | plot_id),
  data = merged_data
)
r1_mod <- glmer.nb(
  aphidius_colemani ~ augmentation * margin_distance + (1|time_point),
  data = merged_data
)
# RQ2 visualisations ----
# RQ3: wasp -  aphid relationship ----
# incrop aphids response
r3_mod_crop <- glmer.nb(
  aphidius_colemani ~ sum_infested * margin_distance + augmentation + time_point +
    (1 | plot_id),
  data = merged_data
)

# RQ3: wasp - aphid relationship ----
# yst aphid response
r3_mod_yst <- glmer.nb(
  aphidius_colemani ~ yst_aphids * margin_distance + augmentation + time_point +
    (1 | plot_id),
  data = merged_data
)

# RQ3 visualisations ----
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

# Do the same for mid-point cover values

wildflower_mc <- wfs_data %>%
  select(`common name`, Species, ends_with("m_c_r")) %>%
  pivot_longer(cols = -c(`common name`, Species),
    names_to = "quadrat",
    values_to = "midpoint_cover"
    ) %>%
  mutate(quadrat = gsub("_m_c_r", "", quadrat))


# merge the two data sets
wfs_long <- wildflower_bb %>%
    left_join(wildflower_mc, by = c("common name", "Species", "quadrat")) %>%
  arrange(quadrat)

view(wfs_long)

# Step 1: Calculate frequency and total cover per species
species_summary <- wfs_long %>%
  group_by(Species) %>%
  summarise(
    plots_of_occurrence = n_distinct(quadrat),
    total_cover = sum(midpoint_cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent_frequency = (plots_of_occurrence / max(plots_of_occurrence)) * 100,
    relative_frequency = (percent_frequency / sum(percent_frequency)) * 100,
    relative_cover = (total_cover / sum(total_cover)) * 100,
    importance_value = relative_frequency + relative_cover
  )


# Sort by importance value
species_summary <- species_summary %>%
  arrange(desc(importance_value))

# View results
View(species_summary)
# Importance value: Fagopyrum esculentum = 74.74453
# Importance value: Phacelia tanacetifolia = 62.57908
# Importance value: Trifolium incarnatum = 21.62633
# Importance value: Trifolium pratense = 20.69151
# Importance value: Myosotis sylvatica = 20.35856

# Average Braun-Blanquet (midpoint cover) per species across all 10 WFS quadrats
wfs_sum_cover <- wfs_long %>%
  group_by(Species) %>% # i think this might be incorrect but the index value looks much smaller without it
  summarise(
    sum_midpoint_cover = sum(midpoint_cover, na.rm = TRUE),
    .groups = "drop"
  )

# in-plot wildflower data wrangling ----

inplot_long <- plot_wf_data %>%
  mutate(across(starts_with("P"), as.character)) %>%
  pivot_longer(cols = starts_with("P"),
    names_to = "plot_position",
    values_to = "braun_blanquet") %>%
  separate(col = plot_position,
    into = c("plot_id", "column"),  # Quoted column names
    sep = "C") %>%                        # Get rid of text from  "C" (e.g., "P1C2" → "P1")
  mutate(
    plot_id = as.numeric(gsub("P", "", plot_id)),  # Convert "P1" → 1
    column = paste0("C", column)                # Restore quadrat label (e.g., "C2")
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
    ))

view(inplot_long)

quadrat_cover <- inplot_long %>%
  group_by(plot_id, column) %>%  # Group by plot and quadrat (column)
  summarise(total_quadrat_cover = sum(midpoint_cover, na.rm = TRUE), .groups = "drop")

view(quadrat_cover)
plot_cover <- quadrat_cover %>%
  group_by(plot_id) %>%
  summarise(mean_wildflower_cover_pct = mean(total_quadrat_cover, na.rm = TRUE))

view(plot_cover)

# add to merged data 
merged_data2 <- inner_join(
  merged_data,
  plot_cover,
  by = "plot_id")
view(merged_data2)

# model iterations in lme4
mod1 <- glmer.nb(
  mean_colemani ~ aphid_mean_infested * mean_wildflower_cover_pct + (1|plot_id),
  data = merged_data2)
summary(mod1)
plot(mod1)
mod1 <- glmer.nb(
  mean_colemani ~ aphid_mean_infested * mean_wildflower_cover_pct + (1|plot_id),
  data = merged_data2)



# lots of unneccessary inplot data manipulations ----
unique(plot_wf_data$Species) 
library(dplyr)

# Calculate global IVs for all species
global_iv <- inplot_long %>%
  group_by(Species) %>%
  summarise(
    plots_of_occurrence = n_distinct(plot_id),
    total_cover = sum(midpoint_cover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    relative_frequency = (plots_of_occurrence / 16) * 100,
    relative_cover = (total_cover / sum(total_cover)) * 100,
    global_iv = relative_frequency + relative_cover
  )
view(global_iv)
# Assign global IVs to each plot and calculate plot-level score
plot_scores <- inplot_long %>%
  left_join(global_iv %>% select(Species, global_iv), by = "Species") %>%
  group_by(plot_id) %>%
  summarise(
    plot_iv = sum(global_iv, na.rm = TRUE),
    total_cover = sum(midpoint_cover, na.rm = TRUE)
  )




# View plot-level scores
View(plot_scores)


# Average Braun-Blanquet (midpoint cover) per species across 5 quadrats per plot
inplot_avg <- plot_wf_data %>%
  group_by(plot_id, Species) %>%  # Group by plot and species
  summarise(
    avg_midpoint_cover = mean(midpoint_cover, na.rm = TRUE),
    .groups = "drop"
  )
view(inplot_avg)
# Calculate Shannon diversity per plot UNNECCESSARY STUFF!!!!!!!!!!---- 
inplot_diversity <- inplot_avg %>%
  group_by(plot_id) %>%
  summarise(
    shannon_div = diversity(avg_midpoint_cover, index = "shannon")
  )
view(inplot_diversity)










view(plot_wf_data)
str(plot_wf_data)

plot_metrics <- plot_wf_data %>%
  group_by(plot_id) %>%
  summarise(
    total_cover = sum(mc_scaled, na.rm = TRUE),
    shannon_div = vegan::diversity(mc_scaled, index = "shannon")
  )


view(plot_metrics)

wfs_pooled <- wfs_long %>%
  group_by(Species) %>%
  summarise(
    total_midpoint_cover = sum(midpoint_cover, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Shannon diversity for the entire WFS (pooled data) (still giving small value)
wfs_diversity_pooled <- diversity(wfs_pooled$total_midpoint_cover, index = "shannon")
view(wfs_diversity_pooled)
t.test(inplot_diversity$shannon_div, mu = wfs_diversity_pooled)
