
#### Load libraries

library(tidyverse)


#### Input files

diversity_file <- "data/diversity_all.csv"
reads_file <- "data/All_Samples_ReadCount_sixtyone.csv"

diversity_all <- read.csv(diversity_file, check.names = FALSE)
reads_count <- read.csv(reads_file, check.names = FALSE)


#### Sample classification

MA_present_samples <- c(
  "PAJ", "PBS", "PBV", "PDA", "PDE", "PDH", "PDN", "PDU",
  "PGA", "PHK", "PHL", "PIK", "PMA", "PMP", "PNN", "PPA",
  "PPP", "PRO", "PSE", "DPQ", "DNL", "DBH", "DNZ", "DJP",
  "DJW", "DSL", "DAP", "DLG", "DSP", "DPN", "DSD", "DDW",
  "DHK", "DHS", "DRH", "DRS", "DRG"
)

diversity_all <- diversity_all %>%
  mutate(
    Microcystis = if_else(
      Sample %in% MA_present_samples,
      "Present", "Absent"
    ),
    Region = if_else(
      grepl("^P", Sample),
      "NCR", "Delhi"
    )
  )


#### Prepare sequencing depth data

reads_depth <- reads_count %>%
  select(
    Sample = Sample_Code,
    processed_reads_combined
  )


#### NCR sequencing depth comparison between MA and NMA

ncr_reads <- diversity_all %>%
  filter(Region == "NCR") %>%
  select(Sample, Microcystis) %>%
  left_join(reads_depth, by = "Sample")

wilcox.test(
  processed_reads_combined ~ Microcystis,
  data = ncr_reads,
  exact = FALSE
)

ncr_reads %>%
  group_by(Microcystis) %>%
  summarise(
    n = n(),
    Median_reads = median(processed_reads_combined),
    Minimum_reads = min(processed_reads_combined),
    Maximum_reads = max(processed_reads_combined),
    .groups = "drop"
  )


#### Correlation between sequencing depth and species richness

reads_richness <- diversity_all %>%
  select(Sample, Richness) %>%
  left_join(reads_depth, by = "Sample")

cor.test(
  reads_richness$processed_reads_combined,
  reads_richness$Richness,
  method = "spearman",
  exact = FALSE
)
