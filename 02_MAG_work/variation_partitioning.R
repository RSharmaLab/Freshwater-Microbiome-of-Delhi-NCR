#### Load libraries

library(tidyverse)
library(readxl)
library(variancePartition)


#### Import data

magabund <- read_csv("data/drepmag_abun_with_taxonomy.csv")
chemistry <- read_excel("data/Hydrochemical Analysis_Merged.xlsx")

magabund <- magabund %>%
  select(-DHK)


#### Define samples and environmental variables

sample_cols <- c(
  "DPQ", "DNL", "DVK", "DJR", "DHS", "DBH", "DNZ", "DSV",
  "DJP", "DJW", "DSL", "DAP", "DLG", "DSP", "DPN", "DSD",
  "DDW", "DRH", "DRS", "DRG",
  "PAJ", "PBA", "PBI", "PBK", "PBS", "PBV", "PCH", "PDA",
  "PDE", "PDH", "PDI", "PDJ", "PDN", "PDO", "PDR", "PDU",
  "PDV", "PDW", "PDY", "PGA", "PHK", "PHL", "PIK", "PJA",
  "PJJ", "PJM", "PKV", "PMA", "PMB", "PMP", "PNA", "PNK",
  "PNN", "PPA", "PPP", "PRO", "PSA", "PSE", "PSU", "PVL"
)

selected_vars <- c(
  "pH", "BOD", "Total_Hardness", "Nitrate", "Nitrite",
  "Ammoniacal_Nitrogen", "Phosphate", "Copper", "Iron",
  "Zinc", "Temperature", "Chromium"
)

chem_reduced <- chemistry %>%
  select(Sample, all_of(selected_vars))


#### Prepare MAG abundance matrix

mag_abund <- magabund %>%
  select(Genome, all_of(sample_cols))

mag_long <- mag_abund %>%
  pivot_longer(
    -Genome,
    names_to = "Sample",
    values_to = "Abundance"
  )

mag_wide <- mag_long %>%
  pivot_wider(
    names_from = Genome,
    values_from = Abundance
  )


#### Merge abundance and chemistry data

dat_mag <- chem_reduced %>%
  inner_join(mag_wide, by = "Sample")

mag_mat <- dat_mag %>%
  select(-Sample, -all_of(selected_vars))

chem_design <- dat_mag %>%
  select(all_of(selected_vars))

mag_mat_log <- t(log2(as.matrix(mag_mat) + 1e-6))

colnames(mag_mat_log) <- dat_mag$Sample
rownames(chem_design) <- dat_mag$Sample


#### Variance partitioning

form <- ~ pH + BOD + Total_Hardness + Nitrate + Nitrite +
  Ammoniacal_Nitrogen + Phosphate + Copper + Iron + Zinc +
  Temperature + Chromium

varPart <- fitExtractVarPartModel(
  mag_mat_log,
  form,
  chem_design
)

errors <- attr(varPart, "errors")


#### SGB prevalence filtering

mag_meta <- magabund %>%
  select(Genome, Phylum)

sgb_prevalence <- data.frame(
  Genome = colnames(mag_mat),
  N_waterbodies = colSums(as.matrix(mag_mat) > 0, na.rm = TRUE)
) %>%
  left_join(
    mag_meta %>% distinct(Genome, Phylum),
    by = "Genome"
  ) %>%
  mutate(Retained = N_waterbodies >= 10)

sgb_counts_phylum <- sgb_prevalence %>%
  group_by(Phylum) %>%
  summarise(
    SGB_before_filter = n(),
    SGB_after_filter = sum(Retained),
    .groups = "drop"
  ) %>%
  arrange(desc(SGB_before_filter))

sgb_counts_phylum
colSums(sgb_counts_phylum[, c("SGB_before_filter", "SGB_after_filter")])


#### Retained SGBs with successful variance partitioning

sgb_keep <- sgb_prevalence %>%
  filter(Retained) %>%
  pull(Genome)

sgb_keep_varpart <- intersect(sgb_keep, rownames(varPart))

cat("SGBs passing prevalence filter:", length(sgb_keep), "\n")
cat("SGBs with successful variancePartition results:",
    length(sgb_keep_varpart), "\n")
cat("SGBs without variancePartition results:",
    length(setdiff(sgb_keep, rownames(varPart))), "\n")


#### Relative environmental contributions per SGB

varpart_retained <- as.data.frame(
  varPart[sgb_keep_varpart, , drop = FALSE]
) %>%
  rownames_to_column("Genome")

variable_cols <- setdiff(
  colnames(varpart_retained),
  c("Genome", "Residuals")
)

varpart_relative <- varpart_retained %>%
  mutate(
    Variable_total = rowSums(across(all_of(variable_cols)))
  ) %>%
  mutate(
    across(
      all_of(variable_cols),
      ~ (.x / Variable_total) * 100
    )
  ) %>%
  select(Genome, all_of(variable_cols)) %>%
  left_join(
    mag_meta %>% distinct(Genome, Phylum),
    by = "Genome"
  )


#### Mean relative contributions by phylum

varpart_phylum <- varpart_relative %>%
  group_by(Phylum) %>%
  summarise(
    N_SGB = n(),
    across(
      all_of(variable_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(N_SGB >= 3) %>%
  arrange(desc(N_SGB))


#### R2 and adjusted R2 per SGB

lm_results <- lapply(sgb_keep, function(id) {
  
  lm_data <- data.frame(
    Abundance = as.numeric(mag_mat_log[id, ]),
    chem_design
  )
  
  fit <- lm(Abundance ~ ., data = lm_data)
  model_summary <- summary(fit)
  
  data.frame(
    Genome = id,
    R2 = model_summary$r.squared,
    Adjusted_R2 = model_summary$adj.r.squared
  )
  
}) %>%
  bind_rows()


#### Mean R2 and adjusted R2 by phylum

lm_phylum <- lm_results %>%
  left_join(
    mag_meta %>% distinct(Genome, Phylum),
    by = "Genome"
  ) %>%
  group_by(Phylum) %>%
  summarise(
    N_SGB = n(),
    Mean_R2 = mean(R2, na.rm = TRUE),
    Mean_Adjusted_R2 = mean(Adjusted_R2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(N_SGB >= 3) %>%
  arrange(desc(N_SGB))


#### Combine phylum-level results

phylum_results <- lm_phylum %>%
  left_join(varpart_phylum, by = c("Phylum", "N_SGB")) %>%
  select(
    Phylum, N_SGB,
    Mean_R2, Mean_Adjusted_R2,
    all_of(variable_cols)
  ) %>%
  arrange(desc(N_SGB))


#### Dominant environmental variable per phylum

dominant_variable <- phylum_results %>%
  rowwise() %>%
  mutate(
    Dominant_variable = variable_cols[
      which.max(c_across(all_of(variable_cols)))
    ],
    Dominant_contribution = max(c_across(all_of(variable_cols)))
  ) %>%
  ungroup() %>%
  select(
    Phylum, N_SGB,
    Dominant_variable, Dominant_contribution
  ) %>%
  arrange(desc(N_SGB))


#### Export results

dir.create("results", showWarnings = FALSE)

write.csv(
  phylum_results,
  "results/Phylum_level_variance_partitioning_results.csv",
  row.names = FALSE
)

write.csv(
  dominant_variable,
  "results/Dominant_environmental_variable_per_phylum.csv",
  row.names = FALSE
)


#### Final summaries

cat("Total SGBs analysed:", nrow(lm_results), "\n")
cat("SGBs retained after prevalence filtering:", length(sgb_keep), "\n")
cat("Phyla included in final analysis:", nrow(phylum_results), "\n")