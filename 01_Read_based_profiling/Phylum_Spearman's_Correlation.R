
#### Load required packages

library(dplyr)
library(tidyr)
library(readr)
library(tibble)

#### Input and output paths

input_file <- "Combined_Rural_Urban_Phylum_Abundance.csv"
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#### Import data

phylum <- read_csv(input_file, show_col_types = FALSE)

#### Phylum abundance comparison between Delhi and NCR

samples <- colnames(phylum)[-1]

region_info <- tibble(
  Sample = samples,
  Region = if_else(grepl("^P", samples), "NCR", "Delhi")
)

phylum_long <- phylum %>%
  pivot_longer(
    cols = -Phylum,
    names_to = "Sample",
    values_to = "Relative_abundance"
  ) %>%
  left_join(region_info, by = "Sample")

selected_phyla <- c("Cyanobacteria", "Proteobacteria")

#### Summary statistics

phylum_summary <- phylum_long %>%
  filter(Phylum %in% selected_phyla) %>%
  group_by(Phylum, Region) %>%
  summarise(
    n = sum(!is.na(Relative_abundance)),
    Median = median(Relative_abundance, na.rm = TRUE),
    Minimum = min(Relative_abundance, na.rm = TRUE),
    Maximum = max(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  phylum_summary,
  file.path(output_dir, "Phylum_Abundance_Summary.csv")
)

#### Wilcoxon tests and rank-biserial correlations

phylum_tests <- phylum_long %>%
  filter(Phylum %in% selected_phyla) %>%
  group_by(Phylum) %>%
  group_modify(~ {
    
    dat <- .x %>%
      filter(!is.na(Relative_abundance))
    
    test <- wilcox.test(
      Relative_abundance ~ Region,
      data = dat,
      exact = FALSE
    )
    
    n_NCR <- sum(dat$Region == "NCR")
    n_Delhi <- sum(dat$Region == "Delhi")
    
    rank_sum_NCR <- sum(
      rank(dat$Relative_abundance, ties.method = "average")[
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
  phylum_tests,
  file.path(output_dir, "Phylum_Wilcoxon_Tests.csv")
)

#### Phylum prevalence and Spearman correlations

phylum_prevalence <- phylum %>%
  mutate(
    n_detected = rowSums(
      across(-Phylum, ~ !is.na(.) & . > 0)
    )
  ) %>%
  select(Phylum, n_detected)

write_csv(
  phylum_prevalence,
  file.path(output_dir, "Phylum_Prevalence.csv")
)

#### Retain phyla detected in at least 10 waterbodies

phylum_keep <- phylum_prevalence %>%
  filter(n_detected >= 10) %>%
  pull(Phylum)

phylum_filtered <- phylum %>%
  filter(Phylum %in% phylum_keep)

phylum_mat <- phylum_filtered %>%
  select(-Phylum) %>%
  as.matrix()

rownames(phylum_mat) <- phylum_filtered$Phylum

phylum_mat <- t(phylum_mat)

#### Spearman correlation and BH correction

cor_results <- Hmisc::rcorr(
  phylum_mat,
  type = "spearman"
)

cor_mat <- cor_results$r
p_mat <- cor_results$P

pair_indices <- which(
  upper.tri(p_mat) & !is.na(p_mat),
  arr.ind = TRUE
)

phylum_correlations <- tibble(
  Phylum_1 = rownames(cor_mat)[pair_indices[, 1]],
  Phylum_2 = colnames(cor_mat)[pair_indices[, 2]],
  Rho = cor_mat[pair_indices],
  p_value = p_mat[pair_indices]
) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH")
  )

write_csv(
  phylum_correlations,
  file.path(output_dir, "Phylum_Correlation_All.csv")
)

#### Retain significant correlations

significant_correlations <- phylum_correlations %>%
  filter(
    abs(Rho) >= 0.5,
    p_adjusted <= 0.05
  ) %>%
  arrange(p_adjusted)

write_csv(
  significant_correlations,
  file.path(output_dir, "Phylum_Correlation_Significant.csv")
)
