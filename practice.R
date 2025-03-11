# RMD
# Diss data analysis
# load packages
library(tidyverse)
library(glmmTMB)
library(MuMIn)
#read in csv
yst_data <- read_excel("data/RMD_sticky_trap_data.xlsx")
view(yst_data)

str(yst_data)

yst_data <- yst_data %>% 
  mutate(across(c(sample_id, set_no, n_treatment, crop, margin, con_rel), as.factor))

mean(yst_data$aphidius_colemani)   #  mean
var(yst_data$aphidius_colemani) # variance

# trialling differnet models and assessing AIC
nb_m1 <- glmmTMB(aphidius_colemani~con_rel * margin + (1 | sample_id) + (1 | set_no),
  family = nbinom2(),  # Variance = μ + μ²/θ
  data = yst_data)

nb_m1 <- glmmTMB(aphidius_colemani~con_rel * margin + (1 | sample_id) + (1 | set_no),
                   family = nbinom2(),  # Variance = μ + μ²/θ
                   data = yst_data)
AICc(nb_m1)
