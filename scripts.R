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
library(readxl) # interpret spreadsheets
library(tidyverse) # dplyr, ggplot, etc.
library(lme4) # preliminary models
library(glmmTMB) # models
library(MuMIn)# AICc values
library(DHARMa) # model diagnostics 
library(ggeffects) # model visualisations
library(patchwork)# multi-plot model visualisations
library(vegan) # For diversity indices? - seems heavily biased for this type of research

#read in data sets
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

wasp_aggregated <- wasp_data %>% 
      group_by(plot_id, time_point) %>%
  arrange(plot_id)
#view(wasp_aggregated)

aphid_aggregated <- aphid_data %>%
  group_by(plot_id, time_point) %>%
  summarise(
    #aphid_total_infested = sum(no_infested),  # Total across all time/sample points
    sum_infested = sum(no_infested),  # Mean per sample point
    .groups = "drop")
#view(aphid_aggregated)

# Combine data sets
merged_data <- inner_join(
  aphid_aggregated, 
  wasp_aggregated, 
  by = c("plot_id", "time_point"))
#str(wasp_aggregated)
view(merged_data)
str(merged_data)
# preliminary analysis ----
# assess distributions of response and predictor variables
hist(merged_data$aphidius_colemani)
shapiro.test(merged_data$aphidius_colemani)

hist(merged_data$sum_infested)
shapiro.test(merged_data$sum_infested)

hist(merged_data$aphidius_colemani)
shapiro.test(merged_data$aphidius_colemani)

hist(merged_data$aphidius_colemani)
shapiro.test(merged_data$aphidius_colemani)

# 
check_overdispersion <- function(response_var) {
  var <- merged_data[[response_var]]
  dispersion <- var(var) / mean(var)
  cat(response_var, "dispersion ratio:", dispersion, "\n")
}

check_overdispersion("aphidius_colemani")  # 1.98 > 1.5 so over dispersed 
check_overdispersion("sum_infested") # 7.2 is massively over dispersed

table(merged_data$aphidius_colemani)  
table(merged_data$sum_infested)  
table(merged_data$yst_aphids)  
# preliminary linear regressions to 
lm1 <-lm(aphidius_colemani ~ sum_infested, data = merged_data) # explains no variation -R2
summary(lm1)
lm2 <-lm(aphidius_colemani ~ augmentation, data = merged_data) # explains almost no variation - R2
summary(lm2)
lm3 <-lm(aphidius_colemani ~ margin_distance, data = merged_data) # explains 5% variation R2 = 0.05
summary(lm3)
lm4 <-lm(aphidius_colemani ~ time_point, data = merged_data) # explains almost no variation -R2
summary(lm4)
lm5 <-lm(aphidius_colemani ~ plot_id, data = merged_data) # explains 25% variation R2 = 0.25, p = 0.045 
summary(lm5)
# incrop aphids
lm6 <-lm(sum_infested ~ augmentation, data = merged_data) # explains  no variation 
summary(lm6)
lm7 <-lm(sum_infested ~ margin_distance, data = merged_data) # explains 4%  variation 
summary(lm7)
lm8 <-lm(sum_infested ~ time_point, data = merged_data) # explains 5%  variation 
summary(lm8)
lm9 <-lm(sum_infested ~ plot_id, data = merged_data) # explains almost no variation 
summary(lm9)
# yst aphids
lm10 <-lm(yst_aphids ~ augmentation, data = merged_data) # explains almost no variation 
summary(lm10)
lm11 <-lm(yst_aphids ~ margin_distance, data = merged_data) # explains almost no variation 
summary(lm11)
lm12 <-lm(yst_aphids ~ time_point, data = merged_data) # explains almost no variation 
summary(lm12)
lm13 <-lm(yst_aphids ~ plot_id, data = merged_data) # explains almost no variation 
summary(lm13)

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

# create themes/palettes  ----
  accent = "#4682B4"          # Steel Blue (optional tertiary)

# RQ1: augmentation and wildflowers effect on colemani ----

r1_mod_yst <- glmmTMB(
  aphidius_colemani ~ augmentation * margin_distance + yst_aphids + (1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)
r1_mod_crop <- glmmTMB(
  aphidius_colemani ~ augmentation * margin_distance + sum_infested + (1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)

check_collinearity(r1_mod_yst)  # For models without interactions
simulateResiduals(r1_mod_yst) %>% plot()

testDispersion(r1_mod_yst)
testZeroInflation(r1_mod_yst)
testUniformity(r1_mod_yst)
testOutliers(r1_mod_yst)

testDispersion(r1_mod_crop)
testZeroInflation(r1_mod_crop)
testUniformity(r1_mod_crop)
testOutliers(r1_mod_crop)

summary(r1_mod_yst)  
summary(r1_mod_crop)  
AICc(r1_mod_yst, r1_mod_crop) # 187.13, 195

# RQ1 visualisations ----
# Get predicted marginal means
aug_margin_effects <- ggpredict(
  r1_mod_crop,
  terms = c("augmentation", "margin_distance", "sum_infested"))
aug_margin_effects <- ggpredict(
  r1_mod_yst,
  terms = c("augmentation", "margin_distance", "yst_aphids"))
view(aug_margin_effects)

# Generate predictions for the interaction
pred_interaction <- ggpredict(
  r1_mod_yst,
  terms = c("augmentation", "margin_distance")  # Interaction only
)
view(pred_interaction)
# Interaction plot
# Generate predictions for the interaction
pred_interaction <- ggpredict(
  r1_mod_yst,
  terms = c("augmentation", "margin_distance")
)

# Plot raw data (boxplots) and model predictions (points ± CI)
p_interaction <- ggplot(merged_data,aes(x = augmentation, y = aphidius_colemani, fill = margin_distance)) + # Raw data boxplots
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_point(data = pred_interaction,
    aes(x = x, y = predicted, group = group),
    color = "black",
    size = 3,
    shape = 18,  # Diamond shape for predictions
    position = position_dodge(width = 0.8)
  ) +
  geom_errorbar(
    data = pred_interaction,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group),
    width = 0.1,
    color = "black",
    linewidth = 1,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_brewer(
    palette = "Dark2", 
    labels = c("33 m", "83 m"),
    name = "WFS Distance"
  ) +
  scale_x_discrete(labels = c("No Augmentation", "Augmentation")) +
  labs(
    x = "Augmentation",
    y = "Aphidius colemani Abundance"
  ) +
  theme_classic(base_size
p_interaction
final_plot <- (p_aug + p_dist) + p_interaction +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

final_plot
# single predictor visualisations

pred_aug <- ggpredict(r1_mod_yst, terms = "augmentation")# Augmentation effect
pred_dist <- ggpredict(r1_mod_yst, terms = "margin_distance")# Distance effect

# Plot 1: Augmentation effect
p_aug <- ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#1b9e77",
    alpha = 0.6) +
  geom_point(
    data = pred_aug,
    aes(x = x, y = predicted),
    size = 3,
    color = "#1b9e77"
  ) +
  geom_errorbar(
    data = pred_aug,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "#1b9e77",
    linewidth = 1) +
  labs(x = "Augmentation", y = "Wasp Abundance") +
  scale_x_discrete(labels = c("No Augmentation", "Augmentation")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))
p_aug
# Plot 2: Distance effect
p_dist <- ggplot(merged_data, aes(x = margin_distance, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#DF6D16",
    alpha = 0.6) +
  geom_point(
    data = pred_dist,
    aes(x = x, y = predicted),
    size = 3,
    color = "#DF6D16"
  ) +
  geom_errorbar(
    data = pred_dist,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "#DF6D16",
    linewidth = 1
  ) +
  labs(x = "Distance from Wildflower Strip", y = "Wasp Abundance") +
  scale_x_discrete(labels = c("33 m", "83 m")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))
p_dist
# Combine plots
combined_plot <- p_aug + p_dist +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

combined_plot
# RQ2: Augmentation and wildflower strip effect on aphid population ----
# Model for in-crop aphids
r2_mod_crop <- glmmTMB(
  sum_infested ~ augmentation * margin_distance + aphidius_colemani + (1|time_point)+ (1|plot_id),
  family = nbinom2,
  data = merged_data)

summary(r2_mod_crop)
simulateResiduals(r2_mod_crop) %>% plot()
testDispersion(r2_mod_crop)
testZeroInflation(r2_mod_crop)
testUniformity(r2_mod_crop)
testOutliers(r2_mod_crop)
summary(r2_mod_crop)$sigma  # For glmmTMB models
AICc(r2_mod_crop) # 334.87

# Model for YST aphids
# negative binomial for interaction term
r2_mod_yst <- glmmTMB(
  yst_aphids ~ augmentation * margin_distance  + aphidius_colemani +(1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)
summary(r2_mod_yst)
simulateResiduals(r2_mod_yst) %>% plot()
testDispersion(r2_mod_yst)
testZeroInflation(r2_mod_yst)
testUniformity(r2_mod_yst)
testOutliers(r2_mod_yst)
AICc(r2_mod_yst) # 204.17
# RQ2 visualisations ----
# marginal means
r2_pred_aug_crop <- ggpredict(r2_mod_crop,
  terms = c("augmentation"))
summary(r2_mod_crop)
# single predictor visualisations
pred_dist_crop <- ggpredict(r2_mod_crop, terms = "augmentation")# Augmentation effect
r2_pred_aug_crop <- ggpredict(r2_mod_crop, terms = c("augmentation", "margin_distance"))
r2_pred_aug_yst <- ggpredict(r2_mod_yst, terms = c("augmentation", "margin_distance"))
pred_dist_yst <- ggpredict(r2_mod_yst, terms = "margin_distance")# Distance effect

# interaction effect predictors
aug_margin_effects_crop <- ggpredict(
  r2_mod_crop,
  terms = c("augmentation", "margin_distance", "aphidius_colemani"))
aug_margin_effects_yst <-ggpredict(
  r2_mod_yst,
  terms = c("augmentation", "margin_distance", "aphidius_colemani"))
# Plot 1: Augmentation effect incrop
p_aug <- ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#1b9e77",
    alpha = 0.6) +
  geom_point(
    data = r2_pred_aug_crop,
    aes(x = x, y = predicted),
    size = 3,
    color = "#1b9e77"
  ) +
  geom_errorbar(
    data = r2_pred_aug_crop,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "#1b9e77",
    linewidth = 1) +
  labs(x = "Augmentation ", y = "Aphis fabae Abundance") +
  scale_x_discrete(labels = c("Present", "Absent")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))
p_aug
# Plot 2: Distance effect incrop
p_dist <- ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#DF6D16",
    alpha = 0.6) +
  geom_point(
    data = pred_dist_crop,
    aes(x = x, y = predicted),
    size = 3,
    color = "#DF6D16") +
  geom_errorbar(
    data = pred_dist_crop,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "#DF6D16",
    linewidth = 1) +
  labs(x = "Augmentation", y = "Aphis Fabae Plant Infestations") +
  scale_x_discrete(labels = c("No Augmentation", "Augmentation")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))
p_dist
# Combine plots
combined_plot <- p_aug + p_dist +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

combined_plot
# Plot 3: Augmentation effect yst
p_aug_yst <- ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#1b9e77",
    alpha = 0.6) +
  geom_point(
    data = r2_pred_aug_yst,
    aes(x = x, y = predicted),
    size = 3,
    color = "#1b9e77"
  ) +
  geom_errorbar(
    data = r2_pred_aug_yst,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "#1b9e77",
    linewidth = 1) +
  labs(x = "Augmentation", y = "Wasp Abundance") +
  scale_x_discrete(labels = c("No Augmentation", "Augmentation")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))
p_aug_yst
# Plot 2: Distance effect yst
p_dist_yst <- ggplot(merged_data, aes(x = margin_distance, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#DF6D16",
    alpha = 0.6) +
  geom_point(
    data = pred_dist_yst,
    aes(x = x, y = predicted),
    size = 3,
    color = "#DF6D16"
  ) +
  geom_errorbar(
    data = pred_dist_yst,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "#DF6D16",
    linewidth = 1
  ) +
  labs(x = "Distance from Wildflower Strip (m)", y = "Wasp Abundance") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"))
p_dist
# in-crop aphid abundance influenced by augmentation and wfs distance 
ggplot(r2_pred_crop, aes(x = x, y = predicted, fill = group)) +
  geom_boxplot(
    size = 0.5,
    color = "black",
    alpha = 0.8,
    outlier.shape = 8
  ) +
  labs(
    x = "Augmentation",
    y = "Predicted Aphis fabae Abundance",
    fill = "Distance from Wildflower Strip"
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("33 m", "83 m")
  ) +
  scale_x_discrete(
    labels = c("No Augmentation", "Augmentation")) +  # Rename x-axis labels
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "bottom")
# YST  aphid abundance influenced by augmentation and wfs distance 
ggplot(r2_pred_yst, aes(x = x, y = predicted, fill = group)) +
  geom_boxplot(
    size = 0.5,
    color = "black",
    alpha = 0.8,
    outlier.shape = 8
  ) +
  labs(
    x = "Augmentation",
    y = "Predicted Aphis fabae Abundance",
    fill = "Distance from Wildflower Strip"
  ) +
  scale_fill_brewer(
    palette = "Dark2",
    labels = c("33 m", "83 m")
  ) +
  scale_x_discrete(
    labels = c("No Augmentation", "Augmentation")) +  # Rename x-axis labels
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "bottom")

# RQ3: wasp -  aphid relationship ----
#  incrop aphid data (same models)
r3_mod_incrop <- glmmTMB(
  aphidius_colemani ~  augmentation + margin_distance + sum_infested + (1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)
# yst aphid data (same models)
r3_mod_yst <- glmmTMB(
  aphidius_colemani ~ augmentation * margin_distance + yst_aphids + (1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)
# incrop diagnostics
summary(r3_mod_incrop)
AICc(r3_mod_incrop) # AICC slightly lower WITHOUT interaction
simulateResiduals(r3_mod_incrop) %>% plot()
testDispersion(r3_mod_incrop)
testZeroInflation(r3_mod_incrop)
testUniformity(r3_mod_incrop)
testOutliers(r3_mod_incrop)
# yst diagnostics
summary(r3_mod_yst)
AICc(r3_mod_yst)# AICC slightly higher WITH interaction
simulateResiduals(r3_mod_yst) %>% plot()
testDispersion(r3_mod_yst)
testZeroInflation(r3_mod_yst)
testUniformity(r3_mod_yst)
testOutliers(r3_mod_yst)

# RQ3 visualisations ----
# Generate predictions conditioned on augmentation
pred_rq3 <- ggpredict(
  r3_mod_incrop,
  terms = c("sum_infested [all]", "margin_distance"))  # x-axis and grouping

# Plot 1: YST aphids vs. wasps
p1 <- ggplot(merged_data, aes(yst_aphids, aphidius_colemani)) +
  geom_point(color = "#1b9e77", alpha = 0.7) +
  geom_smooth(method = "loess", color = "#d95f02", se = TRUE) +
  labs(x = "YST Aphid Abundance", y = "A. colemani Abundance") +
  theme_classic()

# Plot 2: In-crop aphids vs. wasps
p2 <- ggplot(merged_data, aes(sum_infested, aphidius_colemani)) +
  geom_point(color = "#1b9e77", alpha = 0.7) +
  geom_smooth(method = "loess", color = "#d95f02", se = TRUE) +
  labs(x = "In-Crop Aphid Abundance", y = "") +
  theme_classic()
p1
p2
# Combine plots
final_rq3_plot <- p1 + p2 + 
  plot_annotation(tag_levels = "A")
final_rq3_plot

spearman

# 2: colemani - incrop aphid visualisation - AUGMENTATION EFFECT
pred2_rq3 <- ggpredict(
  r3_mod_incrop,
  terms = c("sum_infested [all]", "augmentation"))


# 3: colemani - yst aphid visualisation - MARGIN DISTANCE EFFECT
pred3_rq3 <- ggpredict(
  r3_mod_yst,
  terms = c("yst_aphids [all]", "margin_distance"))


# 4: colemani - yst aphid visualisation - AUGMENTATION EFFECT
pred3_rq3 <- ggpredict(
  r3_mod_yst,
  terms = c("yst_aphids [all]", "augmentation"))


# random bullshit go ----
# Generate predictions for each plot
pred_yst_aug <- ggpredict(r3_mod_yst, terms = c("yst_aphids", "augmentation"))
pred_yst_dist <- ggpredict(r3_mod_yst, terms = c("yst_aphids", "margin_distance"))
pred_incrop_aug <- ggpredict(r3_mod_incrop, terms = c("sum_infested", "augmentation"))
pred_incrop_dist <- ggpredict(r3_mod_incrop, terms = c("sum_infested", "margin_distance"))
# 1: colemani - incrop aphid visualisation - MARGIN DISTANCE EFFECT
p1 <- ggplot(merged_data, aes(x = sum_infested, y = aphidius_colemani)) +
  geom_point(aes(color = margin_distance), alpha = 0.5, size = 2) +  
  geom_ribbon(
    data = pred_incrop_dist,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2) +
  geom_line(
    data = pred_incrop_dist,
    aes(x = x, y = predicted, color = group),
    linewidth = 1) +
  scale_color_manual(
    name = "Distance from WFS",
    values = c("#DF6D16", "#1b9e77"),  # Orange = 33m, Teal = 83m
    labels = c("33m", "83m")) +
  scale_fill_manual(
    name = "Distance from WFS",
    values = c("#DF6D16", "#1b9e77"),
    labels = c("33m", "83m")) +
  labs(
    x = "Aphis fabae crop infestations per plot",
    y = "Aphidius colemani Abundance") +
  theme_classic() +
  theme(legend.position = "bottom")
p1
# 2: colemani - incrop aphid visualisation - AUGMENTATION EFFECT ----
p2 <- ggplot(merged_data, aes(x = sum_infested, y = aphidius_colemani)) +
  geom_point(aes(color = augmentation), alpha = 0.5) +   # Raw data points
  geom_ribbon(                                           # Model predictions with 95% CI
    data = pred2_rq3,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2) +
  geom_line(
    data = pred2_rq3,
    aes(x = x, y = predicted, color = group),
    linewidth = 1
  ) +
  scale_color_manual(values = c("#DF6D16", "#1b9e77")) +
  scale_fill_manual(values = c("#DF6D16", "#1b9e77")) +
  labs(
    x = "Incrop Aphis fabae Abundance",
    y = "Aphidius colemani Abundance",
    color = "augmentation",
    fill = "augmentation") +
  theme_classic() +
  theme(legend.position = "bottom")
p2
# 3: colemani - yst aphid visualisation - MARGIN DISTANCE EFFECT
p3 <- ggplot(merged_data, aes(x = yst_aphids, y = aphidius_colemani)) +
  geom_point(aes(color = margin_distance), alpha = 0.5) +   # Raw data points
  geom_ribbon(                                              # Model predictions with 95% CI
    data = pred3_rq3,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2) +
  geom_line(
    data = pred3_rq3,
    aes(x = x, y = predicted, color = group),
    linewidth = 1
  ) +
  scale_color_manual(values = c("#DF6D16", "#1b9e77")) +
  scale_fill_manual(values = c("#DF6D16", "#1b9e77")) +
  labs(
    x = "YST Aphis fabae Abundance",
    y = "Aphidius colemani Abundance",
    color = "margin_distance",
    fill = "margin_distance"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
p3
# 4: colemani - yst aphid visualisation - AUGMENTATION EFFECT
p4 <- ggplot(merged_data, aes(x = yst_aphids, y = aphidius_colemani)) +
  geom_point(aes(color = augmentation), alpha = 0.5) +   # Raw data points
  geom_ribbon(                                              # Model predictions with 95% CI
    data = pred3_rq3,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2) +
  geom_line(
    data = pred3_rq3,
    aes(x = x, y = predicted, color = group),
    linewidth = 1
  ) +
  scale_color_manual(values = c("#DF6D16", "#1b9e77")) +
  scale_fill_manual(values = c("#DF6D16", "#1b9e77")) +
  labs(
    x = "YST Aphis fabae Abundance",
    y = "Aphidius colemani Abundance",
    color = "augmentation",
    fill = "augmentation"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
p4
multi_plot <- (p1 + p2) / (p3 + p4) + 
  plot_annotation(
    title = "Relationship Between Aphid and Parasitoid Wasp Abundance",
    subtitle = "Stratified by Augmentation and Wildflower Strip (WFS) Distance",
    tag_levels = "A"
  ) &
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
multi_plot
# wildflower strip data wrangling ----

# Convert columns to character
wfs_data <- wfs_data %>%
  mutate(across(starts_with("Q") & !ends_with("m_c_r"), as.character)) %>%
  mutate(across())

# create new column of braun-Blanquet scores
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

# calculating wfs species importance values (no longer useful for results)----
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
    importance_value = relative_frequency + relative_cover) %>%
  arrange(desc(importance_value))

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

# plot wheeling wildflower data wrangling ----

inplot_long <- plot_wf_data %>%
  mutate(across(starts_with("P"), as.character)) %>%
  pivot_longer(cols = starts_with("P"),
    names_to = "plot_position",
    values_to = "braun_blanquet") %>%
  separate(col = plot_position,
    into = c("plot_id", "column"), 
    sep = "C") %>%                        # Get rid of text

  mutate(
    midpoint_cover = case_when(
      braun_blanquet == "+" ~ 0.1,
      braun_blanquet == "1" ~ 2.5,
      braun_blanquet == "2" ~ 15,
      braun_blanquet == "3" ~ 37.5,
      braun_blanquet == "4" ~ 62.5,
      braun_blanquet == "5" ~ 87.5,
      TRUE ~ 0))
inplot_long <- inplot_long %>%
  mutate(across(plot_id, as.factor))
  
view(inplot_long)
str(inplot_long)
quadrat_cover <- inplot_long %>%
  group_by(plot_id, column) %>%  # Group by plot and quadrat (column)
  summarise(total_quadrat_cover = sum(midpoint_cover, na.rm = TRUE), .groups = "drop")

view(quadrat_cover)
plot_cover <- quadrat_cover %>%
  group_by(plot_id) %>%
  summarise(mean_wildflower_cover_pct = mean(total_quadrat_cover, na.rm = TRUE))

view(plot_cover)
str(quadrat_cover)
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
