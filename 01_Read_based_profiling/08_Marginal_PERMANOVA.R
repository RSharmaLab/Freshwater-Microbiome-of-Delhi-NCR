#### Load libraries

library(tidyverse)
library(readxl)
library(vegan)


#### Import data

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


#### Prepare explanatory variables

permanova_data <- combinedhydro %>%
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


#### Calculate Bray-Curtis dissimilarity

community_mat <- species_clean %>%
  column_to_rownames("Species") %>%
  as.matrix() %>%
  t()

community_mat <- community_mat[
  permanova_data$Sample,
  ,
  drop = FALSE
]

bray_dist <- vegdist(community_mat, method = "bray")


#### Marginal PERMANOVA

set.seed(123)

permanova_marginal <- adonis2(
  bray_dist ~ Microcystis_abundance + Temperature + Region,
  data = permanova_data,
  permutations = 999,
  by = "margin"
)

permanova_marginal


#### Sequential PERMANOVA: Temperature first

set.seed(123)

permanova_temp_first <- adonis2(
  bray_dist ~ Temperature + Region + Microcystis_abundance,
  data = permanova_data,
  permutations = 999,
  by = "terms"
)

permanova_temp_first


#### Sequential PERMANOVA: Region first

set.seed(123)

permanova_region_first <- adonis2(
  bray_dist ~ Region + Temperature + Microcystis_abundance,
  data = permanova_data,
  permutations = 999,
  by = "terms"
)

permanova_region_first


#### Full marginal PERMANOVA

set.seed(123)

permanova_full_marginal <- adonis2(
  bray_dist ~ Microcystis_abundance + Region +
    Ammoniacal_Nitrogen + Copper + Nitrite + Nitrate +
    Temperature + BOD + Iron + pH +
    Total_Hardness + Chromium + Phosphate + Zinc,
  data = permanova_data,
  permutations = 999,
  by = "margin"
)

permanova_full_marginal