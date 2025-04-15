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
library(glmmTMB) # Generalised Linear Mixed Models
library(MuMIn)# AICc values
library(DHARMa) # model diagnostics 
library(ggeffects) # model visualisations
library(patchwork)# multi-plot figures
library(ggtext) # italics in plots
library(broom.mixed)  # For tidying GLMMs (Appendix tables)
library(purrr)         # For applying functions across models (appendix tables)
library(flextable)     # For creating formatted tables (for appendix)

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
atr(wfs_data)
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

# Combine data sets
merged_data <- inner_join(
  aphid_aggregated, 
  wasp_aggregated, 
  by = c("plot_id", "time_point"))
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
check_overdispersion("yst_aphids") # 1.4 > 1 so is marginally overdispersed

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


# RQ1: augmentation and wildflowers effect on colemani ----

r1_mod_yst <- glmmTMB(
  aphidius_colemani ~ augmentation * margin_distance + yst_aphids + (1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)

r1_mod_crop <- glmmTMB(
  aphidius_colemani ~ augmentation * margin_distance + sum_infested + (1 | plot_id) + (1 | time_point),
  family = nbinom2,
  data = merged_data)
# models suitability ----
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

summary(r1_mod_yst)  # 
summary(r1_mod_crop)  
AICc(r1_mod_yst, r1_mod_crop) # 187.13, 195.72

# RQ1 visualisations ----
# Get predicted marginal means
aug_margin_effects <- ggpredict(
  r1_mod_crop,
  terms = c("augmentation", "margin_distance", "sum_infested"))
aug_margin_effects <- ggpredict(
  r1_mod_yst,
  terms = c("augmentation", "margin_distance", "yst_aphids"))
view(aug_margin_effects)


# Interaction plot
# Generate predictions for the interaction
pred_interaction <- ggpredict(
  r1_mod_yst,
  terms = c("augmentation", "margin_distance")
)

# Plot raw data (boxplots) and model predictions (points ± CI)
p_interaction <- ggplot(
  merged_data,
  aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    aes(fill = margin_distance),  # Moved fill here
    alpha = 0.8,
    position = position_dodge(width = 0.8)) +
  # Model predictions
  geom_point(
    data = pred_interaction,
    aes(x = x, y = predicted, group = group),
    color = "black",
    size = 3,
    shape = 18,
    position = position_dodge(width = 0.8)) +
  geom_errorbar(
    data = pred_interaction,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group),
    width = 0.1,
    color = "black",
    linewidth = 1,
    position = position_dodge(width = 0.8)) +
  scale_fill_brewer(
    palette = "Dark2", 
    labels = c("33", "83"),
    name = "Distance from Wildflower Strip (m)") +
  scale_x_discrete(
    labels = c("Non-augmented", "Augmented")) +
  labs(
    x = "Augmentation",
    y = "*Aphidius colemani* Abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom", 
        axis.title = element_text(face = "bold"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())

p_interaction
# single predictor visualisations
pred_aug <- ggpredict(r1_mod_yst, terms = "augmentation")# Augmentation effect
pred_dist <- ggpredict(r1_mod_yst, terms = "margin_distance")# Distance effect

# Plot 1: Augmentation effect
p_aug <- ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.7) +
  geom_point(
    data = pred_aug,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = pred_aug,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1) +
  labs(x = "Augmentation", y = "*Aphidius colemani* Abundance") +
  scale_x_discrete(labels = c("Non-augmented", "Augmented")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),axis.title.y = element_markdown(),
    axis.title.x = element_markdown())
p_aug
# Plot 2: Distance effect
p_dist <- ggplot(merged_data, aes(x = margin_distance, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.6) +
  geom_point(
    data = pred_dist,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = pred_dist,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1
  ) +
  labs(x = "Distance from Wildflower Strip (m)", y = "*Aphidius colemani* Abundance") +
  scale_x_discrete(labels = c("33", "83")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),axis.title.y = element_markdown(),
    axis.title.x = element_markdown())
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
summary(r3_mod_yst)
summary(r2_mod_crop)
simulateResiduals(r2_mod_crop) %>% plot()
testDispersion(r2_mod_crop)
testZeroInflation(r2_mod_crop)
testUniformity(r2_mod_crop)
testOutliers(r2_mod_crop)
summary(r2_mod_crop)$sigma  # For glmmTMB models
AICc(r2_mod_crop, r2_mod_yst) # 336.34

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
# RQ2 visualisations - Interactions ----
# interaction effect predictors for incrop and yst aphid data
pred_interaction_crop <- ggpredict(
  r2_mod_crop,
  terms = c("augmentation", "margin_distance"))
pred_interaction_yst <-ggpredict(
  r2_mod_yst,
  terms = c("augmentation", "margin_distance"))
view(pred_interaction_yst)
#Plot raw data (boxplots ) and model predictions (geom_point)
p_interaction_crop <- ggplot(
  merged_data,
  aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    aes(fill = margin_distance),  # Moved fill here
    alpha = 0.8,
    position = position_dodge(width = 0.8)) +
  # Model predictions (no fill needed)
  geom_point(
    data = pred_interaction_crop,
    aes(x = x, y = predicted, group = group),
    color = "black",
    size = 3,
    shape = 18,
    position = position_dodge(width = 0.8)) +
  geom_errorbar(
    data = pred_interaction_crop,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group),
    width = 0.1,
    color = "black",
    linewidth = 1,
    position = position_dodge(width = 0.8)) +
  scale_fill_brewer(
    palette = "Dark2", 
    labels = c("33", "83"),
    name = "Distance from Wildflower Strip (m)") +
  scale_x_discrete(
    labels = c("Non-augmented", "Augmented")) +
  labs(
    x = "Augmentation",
    y = "*Aphis fabae* In Crop Infestations") +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom", 
        axis.title = element_text(face = "bold"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())
p_interaction_crop

p_interaction_yst <- ggplot(
  merged_data,
  aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    aes(fill = margin_distance),  # Moved fill here
    alpha = 0.8,
    position = position_dodge(width = 0.8)) +
  # Model predictions (no fill needed)
  geom_point(
    data = pred_interaction_yst,
    aes(x = x, y = predicted, group = group),
    color = "black",
    size = 3,
    shape = 18,
    position = position_dodge(width = 0.8)) +
  geom_errorbar(
    data = pred_interaction_yst,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group),
    width = 0.1,
    color = "black",
    linewidth = 1,
    position = position_dodge(width = 0.8)) +
  scale_fill_brewer(
    palette = "Dark2", 
    labels = c("33", "83"),
    name = "Distance from Wildflower Strip (m)") +
  scale_x_discrete(
    labels = c("Non-augmented", "Augmented")) +
  labs(
    x = "Augmentation",
    y = "*Aphis fabae* Sticky Trap Abundance") +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none", 
        axis.title = element_text(face = "bold"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())
p_interaction_yst

aphid_int_combined <-  p_interaction_yst / p_interaction_crop +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

aphid_int_combined
# RQ2 visualisations -Single variable plot visualisations ----
# single predictor visualisations
r2_pred_dist_crop <- ggpredict(r2_mod_crop, terms = "margin_distance")# Augmentation effect
r2_pred_aug_crop <- ggpredict(r2_mod_crop, terms = c("augmentation"))
r2_pred_aug_yst <- ggpredict(r2_mod_yst, terms = c("augmentation"))
r2_pred_dist_yst <- ggpredict(r2_mod_yst, terms = "margin_distance")# Distance effect

# Plot 1: Augmentation effect incrop
p_aug <- ggplot(merged_data, aes(x = augmentation, y = aphidius_colemani)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.6) +
  geom_point(
    data = pred_aug,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = pred_aug,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1) +
  labs(x = "Augmentation", y = "*Aphidius colemani* Abundance") +
  scale_x_discrete(labels = c("Non-augmented", "Augmented")) +   
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),axis.title.y = element_markdown(),
    axis.title.x = element_markdown())
p_aug

# Plot 1: Augmentation effect (in-crop)
p_aug_crop <- ggplot(merged_data, aes(x = augmentation, y = sum_infested)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.6
  ) +
  geom_point(
    data = r2_pred_aug_crop,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = r2_pred_aug_crop,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1
  ) +
  labs(
    x = "Augmentation",
    y = "*Aphis fabae* Abundance<br>(In-Crop Plants)"
  ) +
  scale_x_discrete(labels = c("Non-augmented", "Augmented")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )
p_aug_crop
# Plot 2: Distance effect incrop
p_dist_crop <- ggplot(merged_data, aes(x = margin_distance, y = sum_infested)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.6
  ) +
  geom_point(
    data = r2_pred_dist_crop,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = r2_pred_dist_crop,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1
  ) +
  labs(
    x = "Distance from Wildflower Strip (m)",
    y = ""
  ) +
  scale_x_discrete(labels = c("33", "83")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )
p_dist_crop

# Combine in-crop plots
combined_crop <- p_aug_crop + p_dist_crop +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))
combined_crop

# Plot 1: Augmentation effect (YST)
p_aug_yst <- ggplot(merged_data, aes(x = augmentation, y = yst_aphids)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.6
  ) +
  geom_point(
    data = r2_pred_aug_yst,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = r2_pred_aug_yst,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1
  ) +
  labs(
    x = "Augmentation",
    y = "*Aphis fabae* Abundance<br>(Yellow Sticky Traps)"
  ) +
  scale_x_discrete(labels = c("Non-augmented", "Augmented")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )
p_aug_yst
# Plot 2: Distance effect (YST)
p_dist_yst <- ggplot(merged_data, aes(x = margin_distance, y = yst_aphids)) +
  geom_boxplot(
    width = 0.6,
    fill = "#CB4154",
    alpha = 0.6
  ) +
  geom_point(
    data = r2_pred_dist_yst,
    aes(x = x, y = predicted),
    size = 3,
    color = "black"
  ) +
  geom_errorbar(
    data = r2_pred_dist_yst,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    width = 0.1,
    color = "black",
    linewidth = 1
  ) +
  labs(
    x = "Distance from Wildflower Strip (m)",
    y = ""
  ) +
  scale_x_discrete(labels = c("33", "83")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    axis.title.x = element_text(face = "bold"))
p_dist_yst
# Combine YST plots
combined_yst <- p_aug_yst + p_dist_yst +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))
combined_yst
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
    labels = c("Non-augmented", "Augmentation")) +  # Rename x-axis labels
  
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
    labels = c("Non-augmented", "Augmentation")) +  # Rename x-axis labels
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "bottom")

# RQ3: wasp -  aphid relationship MODELS----
#  incrop aphid data (same models)
# For A. colemani ~ aphid abundance, testing interactions

# For A. fabae ~ parasitoid abundance, testing interactions
r3_mod_yst1 <- glmmTMB(yst_aphids ~ augmentation * aphidius_colemani + margin_distance * aphidius_colemani +
                         (1 | plot_id) + (1 | time_point),
                       family = nbinom2,
                       data = merged_data)
r3_mod_yst2 <- glmmTMB(aphidius_colemani ~ augmentation * yst_aphids + margin_distance * yst_aphids +
                         (1 | plot_id) + (1 | time_point),
                       family = nbinom2,
                       data = merged_data)
r3_mod_crop1 <- glmmTMB(sum_infested ~ augmentation * aphidius_colemani + margin_distance * aphidius_colemani +
                          (1 | plot_id) + (1 | time_point),
                        family = nbinom2,
                        data = merged_data)
r3_mod_crop2 <- glmmTMB(aphidius_colemani ~ augmentation * sum_infested + margin_distance * sum_infested +
                          (1 | plot_id) + (1 | time_point),
                        family = nbinom2,
                        data = merged_data)

# incrop diagnostics
summary(r3_mod_crop1)
summary(r3_mod_yst1)
summary(r3_mod_crop2)
summary(r3_mod_yst2)
AICc(r3_mod_crop1,r3_mod_yst1,r3_mod_yst2,r3_mod_crop2) 
simulateResiduals(r3_mod_crop1) %>% plot()
simulateResiduals(r3_mod_yst1) %>% plot() # we going with this to visualise
simulateResiduals(r3_mod_crop2) %>% plot() # and this one
simulateResiduals(r3_mod_yst2) %>% plot()

testDispersion(r3_mod_yst1)
testZeroInflation(r3_mod_yst1)
testUniformity(r3_mod_yst1)
testOutliers(r3_mod_yst1)
testDispersion(r3_mod_crop2)
testZeroInflation(r3_mod_crop2)
testUniformity(r3_mod_crop2)
testOutliers(r3_mod_crop2)
# RQ3 visualisations ----
# Predictions for augmentation effect
pred_yst1_aug <- ggpredict(r3_mod_yst1,
                           terms = c("aphidius_colemani", "augmentation"))  # Vary parasitoid abundance and augmentation
# Predictions for WFS proximity effect
pred_yst1_dist <- ggpredict(r3_mod_yst1,
                            terms = c("aphidius_colemani", "margin_distance"))  # Vary parasitoid abundance and WFS distance

# Predictions for augmentation effect
pred_crop2_aug <- ggpredict(r3_mod_crop2,
                            terms = c("sum_infested", "augmentation")) # Vary aphid abundance and augmentation
# Predictions for WFS proximity effect
pred_crop2_dist <- ggpredict(r3_mod_crop2,
                             terms = c("sum_infested", "margin_distance"))  # Vary aphid abundance and WFS distance


# Panel A: Augmentation effect on crop infestations
p_aug_crop <- ggplot(merged_data, aes(x = sum_infested, y = aphidius_colemani)) +
  # Raw data points
  geom_point(aes(color = augmentation), alpha = 0.6, size = 2) +
  # Model predictions with 95% CI
  geom_ribbon(
    data = pred_crop2_aug,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2
  ) +
  geom_line(
    data = pred_crop2_aug,
    aes(x = x, y = predicted, color = group),
    linewidth = 1.2
  ) +
  scale_color_manual(
    name = "Augmentation",
    values = c("#CB4154", "#87CEEB"),  # Red = Non-augmented, Light blue = Augmented
    labels = c("Non-augmented", "Augmented")
  ) +
  scale_fill_manual(
    name = "Augmentation",
    values = c("#CB4154", "#87CEEB"),
    labels = c("Non-augmented", "Augmented")
  ) +
  labs(
    x = "*Aphis fabae* Abundance <br> (Crop Infestations)",
    y = "*Aphidius colemani* Abundance",
    tag = "A"
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_markdown(face = "bold"),
    axis.title.x = element_markdown(face = "bold"),
    legend.position = "none")

# Panel B: Augmentation effect on YST aphids
p_aug_yst <- ggplot(merged_data, aes(x = yst_aphids, y = aphidius_colemani)) +
  geom_point(aes(color = augmentation), alpha = 0.6, size = 2) +
  geom_ribbon(
    data = pred_yst1_aug,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2
  ) +
  geom_line(
    data = pred_yst1_aug,
    aes(x = x, y = predicted, color = group),
    linewidth = 1.2
  ) +
  scale_color_manual(
    name = "Augmentation",
    values = c("#CB4154", "#87CEEB"),
    labels = c("Non-augmented", "Augmented")
  ) +
  scale_fill_manual(
    name = "Augmentation",
    values = c("#CB4154", "#87CEEB"),
    labels = c("Non-augmented", "Augmented")
  ) +
  labs(
    x = "*Aphis fabae* Abundance <br> (Yellow Sticky Traps)",
    y = "",
    tag = "B"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(face = "bold"),
    legend.position = "bottom"
  )

# Combine panels into one figure
aug_figure <- (p_aug_crop | p_aug_yst) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A",theme = theme(legend.position = "bottom"))
aug_figure
#  WFS distance plots:

# A: colemani - incrop aphid visualisation - MARGIN DISTANCE EFFECT
p_dist_crop<- ggplot(merged_data, aes(x = sum_infested, y = aphidius_colemani)) +
  geom_point(aes(color = margin_distance), alpha = 0.5, size = 2) +  
  geom_ribbon(
    data = pred_crop2_dist,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2) +
  geom_line(
    data = pred_crop2_dist,
    aes(x = x, y = predicted, color = group),
    linewidth = 1) +
  scale_color_manual(
    name = "Distance from Wildflower Strip (m)",
    values = c( "#1b9e77","#DF6D16"),  # Orange = 33m, Teal = 83m
    labels = c("33", "83")) +
  scale_fill_manual(
    name = "Distance from Wildflower Strip (m)",
    values = c( "#1b9e77","#DF6D16"),
    labels = c("33", "83")) +
  labs(
    x = "*Aphis fabae* Abundance<br>(Crop Infestations)",
    y = "*Aphidius colemani* Abundance") +
  theme_classic() +
  theme(legend.position = "bottom",axis.title.y = element_markdown(face = "bold"),
        axis.title.x = element_markdown(face = "bold"))
p_dist_crop

# B: colemani - yst aphid visualisation - MARGIN DISTANCE EFFECT
p_dist_yst <- ggplot(merged_data, aes(x = yst_aphids, y = aphidius_colemani)) +
  geom_point(aes(color = margin_distance), alpha = 0.5, size = 2) +   
  geom_ribbon(                                             
    data = pred_yst1_dist,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2) +
  geom_line(
    data = pred_yst1_dist,
    aes(x = x, y = predicted, color = group),
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(0, 30)) +  # Zoom in
  scale_color_manual(values = c("#1b9e77","#DF6D16")) +
  scale_fill_manual(values = c("#1b9e77","#DF6D16")) +
  labs(
    x = "*Aphis fabae* Abundance <br> (Yellow Sticky Traps)",
    y = "",
    color = "margin_distance",
    fill = "margin_distance"
  ) +
  theme_classic() +
  theme(legend.position = "none",axis.title.y = element_markdown(face = "bold"),
        axis.title.x = element_markdown(face = "bold"))
p_dist_yst
# combined WFS distance plots
dist_figure <- (p_dist_crop | p_dist_yst) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A",theme = theme(legend.position = "bottom"))
dist_figure

#RQ4:  plot wheeling wildflower influence on species abundance ----
# wrangling
inplot_long <- plot_wf_data %>%
  mutate(across(starts_with("P"), as.character)) %>%
  pivot_longer(cols = starts_with("P"),
               names_to = "plot_position",
               values_to = "braun_blanquet") %>%
  separate(col = plot_position,
           into = c("plot_id", "column"), 
           sep = "C") %>%                        # Get rid of text and separate into 2 variables
  
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
species_summary <- inplot_long %>%
  group_by(Species) %>%
  summarise(
    Total_Cover_Score = sum(plot_total_cover, na.rm = TRUE),    # Total cover across all plots (sum of plot_total_cover)
    Mean_Cover_Score = round(mean(plot_mean_cover[plot_presence == 1], na.rm = TRUE), 2),    # Mean cover per plot (average of plot_mean_cover, excluding plots where species was Non-augmented)
    Relative_Cover = Total_Cover_Score / sum(Total_Cover_Score),    # Relative cover (proportion of total cover across all species)
    
    .groups = "drop"
  ) %>%
  # Rename columns for clarity
  rename(
    "Species" = Species,
    "Number of Plots" = Number_of_Plots,
    "Total Cover Score" = Total_Cover_Score,
    "Mean Cover Score" = Mean_Cover_Score,
    "Relative Cover" = Relative_Cover)
# get total quadrat cover of all flowers 
quadrat_cover <- inplot_long %>%
  group_by(plot_id, column) %>%  # Group by plot and quadrat (column)
  summarise(total_quadrat_cover = sum(midpoint_cover, na.rm = TRUE), .groups = "drop")

view(quadrat_cover)
# get average per plot
plot_cover <- quadrat_cover %>%
  group_by(plot_id) %>%
  summarise(p_quadrat_cover = (sum(total_quadrat_cover, na.rm = TRUE)/5))

view(plot_cover)
str(plot_cover)
str(merged_data)
# add to merged data 
# Remove "P" from plot_id and convert to factor
plot_cover$plot_id <- gsub("P", "", plot_cover$plot_id)
plot_cover$plot_id <- as.factor(plot_cover$plot_id)

# Merge datasets
merged_data2 <- merge(merged_data, plot_cover, by = "plot_id", all.x = TRUE)
view(merged_data2)

# wasp model iterations ----

wasp_model <- glmmTMB(
  aphidius_colemani ~ p_quadrat_cover + augmentation + margin_distance + 
    (1 | time_point) + (1 | plot_id),
  family = nbinom2,  # Use nbinom2 if overdispersed, else poisson
  data = merged_data2)
summary(wasp_model)
AICc(wasp_model)
simulateResiduals(wasp_model) %>% plot()

# aphid model iterations ----
# For YST aphids
glmm_yst <- glmmTMB(
  yst_aphids ~ p_quadrat_cover + augmentation + margin_distance + (1 | time_point) + (1 | plot_id),
  family = nbinom2,  
  data = merged_data2)

summary(glmm_yst)
AICc(glmm_yst)
simulateResiduals(glmm_yst) %>% plot()

# For sum_infested (in-crop aphids)
glmm_infested <- glmmTMB(
  sum_infested ~ p_quadrat_cover + augmentation + margin_distance + (1 | time_point) + (1 | plot_id),
  family = nbinom2,
  data = merged_data2)
summary(glmm_infested)
AICc(glmm_infested)
simulateResiduals(glmm_infested) %>% plot()
# 1. Wasp Abundance (Aphidius colemani) ----
wasp_pred <- ggpredict(wasp_model, terms = "p_quadrat_cover[0:50]")

plot_wasp <- ggplot() +
  geom_point(data = merged_data2,
             aes(x = p_quadrat_cover, y = aphidius_colemani, color = margin_distance, shape = augmentation),
             size = 2.5, stroke = 0.8, alpha = 0.6 ) +
  geom_line(
    data = wasp_pred,
    aes(x = x, y = predicted),
    color = "navy", linewidth = 1 ) +
  geom_ribbon(
    data = wasp_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "navy", alpha = 0.2 ) +
  labs(
    x = "Wildflower Cover <br>(cm² Index)",
    y = "*Aphidius colemani* Abundance") +
  xlim(0, 50) +
  theme_classic() +
  theme(
    axis.title.y = element_markdown(face = "bold", size = 10),
    axis.title.x = element_markdown(face = "bold", size = 10),
    legend.position = "none" ) +
  # Color scale for margin distance (33m = green, 83m = orange)
  scale_color_manual(
    name = "Distance from Wildflower Strip (m)",
    values = c("33" = "#1b9e77", "83" = "#DF6D16"), 
    labels = c("33", "83") ) +
  # Shape scale for augmentation (present = circle, absent = triangle)
  scale_shape_manual(
    name = "Augmentation",
    values = c("Augmented" = 16, "Non-augmented" = 17),  # 16 = circle, 17 = triangle
    labels = c("Augmented", "Non-augmented"))

plot_wasp
# 2.YST Aphids ----

yst_pred <- ggpredict(glmm_yst, terms = "p_quadrat_cover[0:50]")

plot_yst <- ggplot() +
  geom_point(data = merged_data2,
             aes(x = p_quadrat_cover, y = yst_aphids, color = margin_distance, shape = augmentation),
             size = 2.5, stroke = 0.8, alpha = 0.6 ) +
  geom_line(
    data = yst_pred,
    aes(x = x, y = predicted),
    color = "navy", linewidth = 1 ) +
  geom_ribbon(
    data = yst_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "navy", alpha = 0.2 ) +
  labs(
    x = "Wildflower Cover <br>(cm² Index)",
    y = "*Aphis fabae* Abundance <br> (Yellow Sticky Traps)") +
  xlim(0, 50) +
  theme_classic() +
  theme(
    axis.title.y = element_markdown(face = "bold", size = 10),
    axis.title.x = element_markdown(face = "bold", size = 10),
    legend.position = "none" ) +
  # Color scale for margin distance (33m = green, 83m = orange)
  scale_color_manual(
    name = "Margin Distance",
    values = c("33" = "#1b9e77", "83" = "#DF6D16"), 
    labels = c("33 m", "83 m") ) +
  # Shape scale for augmentation (present = circle, absent = triangle)
  scale_shape_manual(
    name = "Augmentation",
    values = c("Augmented" = 16, "Non-augmented" = 17),  # 16 = circle, 17 = triangle
    labels = c("Augmented", "Non-augmented"))
plot_yst
# 3. Crop Infestations (Aphis fabae in crop) ----
infested_pred <- ggpredict(glmm_infested, terms = "p_quadrat_cover[0:50]")

plot_infested <- ggplot() +
  geom_point(data = merged_data2,
             aes(x = p_quadrat_cover, y = sum_infested, color = margin_distance, shape = augmentation),
             size = 2.5, stroke = 0.8, alpha = 0.6 ) +
  geom_line(
    data = infested_pred,
    aes(x = x, y = predicted),
    color = "navy", linewidth = 1 ) +
  geom_ribbon(
    data = infested_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "navy", alpha = 0.2 ) +
  labs(
    x = "Wildflower Cover <br>(cm² Index)",
    y = "*Aphis fabae* Abundance<br>(Crop Infestations)") +
  xlim(0, 50) +
  theme_classic() +
  theme(
    axis.title.y = element_markdown(face = "bold", size = 10),
    axis.title.x = element_markdown(face = "bold", size = 10),
    legend.position = "none" ) +
  # Color scale for margin distance (33m = green, 83m = orange)
  scale_color_manual(
    name = "Margin Distance",
    values = c("33" = "#1b9e77", "83" = "#DF6D16"), 
    labels = c("33 m", "83 m") ) +
  # Shape scale for augmentation (present = circle, absent = triangle)
  scale_shape_manual(
    name = "Augmentation",
    values = c("present" = 16, "absent" = 17),  # 16 = circle, 17 = triangle
    labels = c("Present", "Absent"))
plot_infested
# 3. Combine Plots with Unified Legends ----
combined_plots <- plot_wasp / (plot_yst + plot_infested) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.tag = element_text(face = "bold", size = 12)
  )

combined_plots

# Final model outputs  for appendicies ----------
models <- list(
  "RQ1: A. colemani Crop (Augmentation * WFS)" = r1_mod_crop,
  "RQ1: A. colemani YST (Augmentation * WFS)" = r1_mod_yst,
  "RQ2: In-Crop Aphids" = r2_mod_crop,
  "RQ2: YST Aphids" = r2_mod_yst,
  "RQ3: YST Aphid-Wasp Interactions" = r3_mod_yst1,  #  YST aphids model
  "RQ3: In-Crop Aphid-Wasp Interactions" = r3_mod_crop2,  #  crop aphids model
  "RQ4: A. colemani & WF correlation" = wasp_model,  #  crop aphids model
  "RQ4: YST A. fabae & WF correlation" = glmm_yst,  #  crop aphids model
  "RQ4: In-Crop A. fabae & WF correlation" = glmm_infested)  #  crop aphids model

# Create summary table
model_table <- map_dfr(models, ~ broom.mixed::tidy(.x, conf.int = TRUE, effects = "fixed"), .id = "Model") %>%
  mutate(across(c(estimate, std.error, statistic, conf.low, conf.high), ~ round(., 3)),
         p.value = round(p.value, 4)) %>%
  select(Model, term, estimate, std.error, statistic, conf.low, conf.high, p.value)

# Get BIC and LogLiklihood statistics
fit_stats <- map_dfr(models, ~ broom.mixed::glance(.x), .id = "Model") %>%
  mutate(across(c( BIC, logLik), ~ round(., 1)))

# Combine tables
full_table <- left_join(model_table, fit_stats, by = "Model")

# View final table
View(full_table)
# Final Model AICcs ----
AICc(r1_mod_yst,r1_mod_crop, r2_mod_crop, r2_mod_yst,r3_mod_crop2,r3_mod_yst1, wasp_model,glmm_yst,glmm_infested )
# Final model diagnostic summary table ----
# RQ1 MODEL 1
testDispersion(r1_mod_yst) #dispersion = 0.97941, p-value = 0.952
testZeroInflation(r1_mod_yst) #ratioObsSim = 1.1966, p-value = 0.576
testUniformity(r1_mod_yst) #D = 0.074571, p-value = 0.9341
testOutliers(r1_mod_yst) #outliers at both margin(s) = 0, observations = 48, p-value = 1

# RQ1 MODEL 2 (subsequent outputs copied directly into excel table)
testDispersion(r1_mod_crop)
testZeroInflation(r1_mod_crop)
testUniformity(r1_mod_crop)
testOutliers(r1_mod_crop)

# RQ2 MODEL 1
testDispersion(r2_mod_crop)
testZeroInflation(r2_mod_crop)
testUniformity(r2_mod_crop)
testOutliers(r2_mod_crop)

# RQ2 MODEL 2
testDispersion(r2_mod_yst)
testZeroInflation(r2_mod_yst)
testUniformity(r2_mod_yst)
testOutliers(r2_mod_yst)

# RQ3 MODEL 1
testDispersion(r3_mod_yst1)
testZeroInflation(r3_mod_yst1)
testUniformity(r3_mod_yst1)
testOutliers(r3_mod_yst1)

# RQ3 MODEL 2
testDispersion(r3_mod_crop2)
testZeroInflation(r3_mod_crop2)
testUniformity(r3_mod_crop2)
testOutliers(r3_mod_crop2)

# RQ4 MODEL 1
testDispersion(wasp_model)
testZeroInflation(wasp_model)
testUniformity(wasp_model)
testOutliers(wasp_model)

# RQ4 MODEL 2
testDispersion(glmm_yst)
testZeroInflation(glmm_yst)
testUniformity(glmm_yst)
testOutliers(glmm_yst)

# RQ4 MODEL 3
testDispersion(glmm_infested)
testZeroInflation(glmm_infested)
testUniformity(glmm_infested)
testOutliers(glmm_infested)
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

# calculating wfs species importance values ----
total_cover_sum <- sum(species_summary$total_cover)

species_summary <- wfs_long %>%
  group_by(Species) %>%
  summarise(
    plots_of_occurrence = n_distinct(quadrat),
    total_cover = sum(midpoint_cover, na.rm = TRUE),
    mean_cover = mean(midpoint_cover, na.rm = T),
    #relative_cover = (total_cover / sum(total_cover)),
    #relative_cover = total_cover / total_cover_sum,
    #relative_cover_percent = scales::percent(relative_cover),  # Optional: format as %
    .groups = "drop") %>%
  rename(
    "Number of Plots of Occurrence" = plots_of_occurrence,  
    "Total Cover Score" = total_cover,  # summed Braun-Blanquet MCR
    "Mean cover (%)" = mean_cover)
#mutate(
#  percent_frequency = (plots_of_occurrence / max(plots_of_occurrence)) * 100,
#  relative_frequency = (percent_frequency / sum(percent_frequency)) * 100,
#  relative_cover = (total_cover / sum(total_cover)) * 100,
#  importance_value = relative_frequency + relative_cover) %>%
#  arrange(desc(importance_value))
# Calculate total cover

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
