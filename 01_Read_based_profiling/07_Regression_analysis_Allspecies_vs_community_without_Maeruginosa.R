#### Load libraries

library(tidyverse)
library(readxl)


#### Import data

diversity_all <- read_csv("data/combined_alpha_diversity_microcystis.csv")
diversity_no_MA <- read_csv("data/diversity_no_MA.csv")
species <- read_csv("data/species_abund_combined.csv")
combinedhydro <- read_excel("data/Hydrochemical Analysis_Merged.xlsx")


#### Extract Microcystis aeruginosa abundance

species_clean <- species %>%
  select(-DHK)

microcystis_abund <- species_clean %>%
  filter(str_detect(Species, "Microcystis aeruginosa")) %>%
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Microcystis_abundance"
  ) %>%
  select(Sample, Microcystis_abundance)


#### Prepare common model data

model_data <- combinedhydro %>%
  filter(Sample != "DHK") %>%
  select(
    Sample, Ammoniacal_Nitrogen, Copper, Nitrite, Nitrate,
    Temperature, BOD, Iron, pH, Total_Hardness, Chromium,
    Phosphate, Zinc
  ) %>%
  left_join(microcystis_abund, by = "Sample") %>%
  mutate(
    Region = factor(
      if_else(str_detect(Sample, "^P"), "NCR", "Delhi"),
      levels = c("Delhi", "NCR")
    )
  )


#### Prepare response datasets

model_data_all <- model_data %>%
  left_join(
    diversity_all %>% select(Sample, Shannon),
    by = "Sample"
  )

model_data_no_MA <- model_data %>%
  left_join(
    diversity_no_MA %>% select(Sample, Shannon),
    by = "Sample"
  )


#### Define regression models

model_A <- Shannon ~ Microcystis_abundance +
  Ammoniacal_Nitrogen + Copper + Nitrite + Nitrate +
  BOD + Iron + pH + Total_Hardness + Chromium +
  Phosphate + Zinc + Region

model_B <- Shannon ~ Microcystis_abundance +
  Ammoniacal_Nitrogen + Copper + Nitrite + Nitrate +
  Temperature + BOD + Iron + pH + Total_Hardness +
  Chromium + Phosphate + Zinc + Region


#### Regression: Original community

lm_shannon_all_A <- lm(model_A, data = model_data_all)
lm_shannon_all_B <- lm(model_B, data = model_data_all)

summary(lm_shannon_all_A)
summary(lm_shannon_all_B)


#### Regression: Community without M. aeruginosa

lm_shannon_no_MA_A <- lm(model_A, data = model_data_no_MA)
lm_shannon_no_MA_B <- lm(model_B, data = model_data_no_MA)

summary(lm_shannon_no_MA_A)
summary(lm_shannon_no_MA_B)


#### Variance inflation factor

car::vif(lm_shannon_all_A)
car::vif(lm_shannon_all_B)

car::vif(lm_shannon_no_MA_A)
car::vif(lm_shannon_no_MA_B)