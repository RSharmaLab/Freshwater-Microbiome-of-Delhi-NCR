
#### Load required packages

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(vegan)

#### Input and output paths

input_file <- "Species_Abundance.csv"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#### Import species abundance data

species <- read_csv(input_file, show_col_types = FALSE)

#### Calculate alpha diversity using all species

species_all <- species %>%
  mutate(
    across(
      -Species,
      ~ (.x / sum(.x, na.rm = TRUE)) * 100
    )
  )

#### Convert to sample-by-species matrix

community_all <- species_all %>%
  column_to_rownames("Species") %>%
  as.matrix() %>%
  t()

#### Calculate Shannon diversity and richness

diversity_all <- data.frame(
  Sample = rownames(community_all),
  Shannon = vegan::diversity(community_all, index = "shannon"),
  Richness = vegan::specnumber(community_all)
)

#### Assign samples to Delhi and NCR

samples <- diversity_all$Sample

region_info <- tibble(
  Sample = samples,
  Region = if_else(grepl("^P", samples), "NCR", "Delhi")
)

alpha_combined <- diversity_all %>%
  left_join(region_info, by = "Sample")

write_csv(
  alpha_combined,
  file.path(output_dir, "Alpha_Diversity_All_Species.csv")
)

#### Wilcoxon tests and rank-biserial effect sizes

alpha_tests <- alpha_combined %>%
  pivot_longer(
    cols = c(Shannon, Richness),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  group_by(Metric) %>%
  group_modify(~ {
    
    dat <- .x %>%
      filter(!is.na(Value))
    
    test <- wilcox.test(
      Value ~ Region,
      data = dat,
      exact = FALSE
    )
    
    n_NCR <- sum(dat$Region == "NCR")
    n_Delhi <- sum(dat$Region == "Delhi")
    
    rank_sum_NCR <- sum(
      rank(dat$Value, ties.method = "average")[
        dat$Region == "NCR"
      ]
    )
    
    U_NCR <- rank_sum_NCR - (n_NCR * (n_NCR + 1)) / 2
    
    r_rb <- (2 * U_NCR / (n_NCR * n_Delhi)) - 1
    
    tibble(
      W = as.numeric(test$statistic),
      p_value = test$p.value,
      Rank_biserial_r = r_rb
    )
  }) %>%
  ungroup()

write_csv(
  alpha_tests,
  file.path(output_dir, "Alpha_Diversity_Wilcoxon_Tests.csv")
)

#### Assign Microcystis presence/absence status based on sample IDs

MA_present_samples <- c(
  "PAJ", "PBS", "PBV", "PDA", "PDE", "PDH", "PDN", "PDU",
  "PGA", "PHK", "PHL", "PIK", "PMA", "PMP", "PNN", "PPA",
  "PPP", "PRO", "PSE", "DPQ", "DNL", "DBH", "DNZ", "DJP",
  "DJW", "DSL", "DAP", "DLG", "DSP", "DPN", "DSD", "DDW",
  "DHK", "DHS", "DRH", "DRS", "DRG"
)

alpha_mapresent <- diversity_all %>%
  mutate(
    Microcystis = if_else(
      Sample %in% MA_present_samples,
      "Present",
      "Absent"
    ),
    Microcystis = factor(
      Microcystis,
      levels = c("Absent", "Present")
    )
  )

write_csv(
  alpha_mapresent,
  file.path(output_dir, "Alpha_Diversity_MA_Presence_Absence.csv")
)

#### Wilcoxon tests and rank-biserial effect sizes

alpha_MA_tests <- alpha_mapresent %>%
  pivot_longer(
    cols = c(Shannon, Richness),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  group_by(Metric) %>%
  group_modify(~ {
    
    dat <- .x %>%
      filter(!is.na(Value))
    
    test <- wilcox.test(
      Value ~ Microcystis,
      data = dat,
      exact = FALSE
    )
    
    n_absent <- sum(dat$Microcystis == "Absent")
    n_present <- sum(dat$Microcystis == "Present")
    
    rank_sum_absent <- sum(
      rank(dat$Value, ties.method = "average")[
        dat$Microcystis == "Absent"
      ]
    )
    
    U_absent <- rank_sum_absent -
      (n_absent * (n_absent + 1)) / 2
    
    r_rb <- 1 - (2 * U_absent / (n_absent * n_present))
    
    tibble(
      W = as.numeric(test$statistic),
      p_value = test$p.value,
      Rank_biserial_r = r_rb
    )
  }) %>%
  ungroup()

write_csv(
  alpha_MA_tests,
  file.path(output_dir, "Alpha_Diversity_MA_Wilcoxon_Tests.csv")
)
