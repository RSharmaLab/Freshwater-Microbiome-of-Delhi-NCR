
#### Load libraries

library(tidyverse)
library(vegan)


#### Input files

species_file <- "data/species_abundance_matrix.csv"
diversity_file <- "data/diversity_all.csv"

species <- read.csv(species_file, check.names = FALSE)
diversity_all <- read.csv(diversity_file, check.names = FALSE)


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


#### NCR-only alpha diversity comparison

ncr_diversity <- diversity_all %>%
  filter(Region == "NCR")


#### Shannon diversity

wilcox.test(
  Shannon ~ Microcystis,
  data = ncr_diversity,
  exact = FALSE
)

n_present <- sum(ncr_diversity$Microcystis == "Present")
n_absent <- sum(ncr_diversity$Microcystis == "Absent")

rank_sum_MA <- sum(
  rank(ncr_diversity$Shannon)[
    ncr_diversity$Microcystis == "Present"
  ]
)

U_MA <- rank_sum_MA - (n_present * (n_present + 1)) / 2

r_rb_shannon <- (2 * U_MA / (n_present * n_absent)) - 1

r_rb_shannon

ncr_diversity %>%
  group_by(Microcystis) %>%
  summarise(
    n = n(),
    Median_Shannon = median(Shannon),
    Minimum = min(Shannon),
    Maximum = max(Shannon),
    .groups = "drop"
  )


#### Species richness

wilcox.test(
  Richness ~ Microcystis,
  data = ncr_diversity,
  exact = FALSE
)

rank_sum_MA <- sum(
  rank(ncr_diversity$Richness)[
    ncr_diversity$Microcystis == "Present"
  ]
)

U_MA <- rank_sum_MA - (n_present * (n_present + 1)) / 2

r_rb_richness <- (2 * U_MA / (n_present * n_absent)) - 1

r_rb_richness

ncr_diversity %>%
  group_by(Microcystis) %>%
  summarise(
    n = n(),
    Median_Richness = median(Richness),
    Minimum = min(Richness),
    Maximum = max(Richness),
    .groups = "drop"
  )


#### NCR microbial community composition

ncr_abundance <- species %>%
  select(starts_with("P")) %>%
  t() %>%
  as.data.frame()

ncr_abundance$Sample <- rownames(ncr_abundance)

ncr_abundance <- ncr_abundance %>%
  left_join(
    diversity_all %>% select(Sample, Microcystis),
    by = "Sample"
  )


#### Bray-Curtis dissimilarity

ncr_bray <- vegdist(
  ncr_abundance %>% select(-Sample, -Microcystis),
  method = "bray"
)


#### PERMANOVA

ncr_permanova <- adonis2(
  ncr_bray ~ Microcystis,
  data = ncr_abundance,
  permutations = 999
)

ncr_permanova

#ALPHA DIVERSITY COMPARISON BETWEEN MA AND NMA SAMPLES IN DELHI

# Subset the alpha diversity dataset to include only Delhi samples.

delhi_diversity <- mashanon %>%
  filter(grepl("^D", Sample))


#SHANNON DIVERSITY

# Compare Shannon diversity between MA and NMA samples
# using the Wilcoxon rank-sum test.

wilcox.test(
  Shannon ~ Microcystis,
  data = delhi_diversity,
  exact = FALSE
)
#p-value = 0.03071

# Group sizes

n_present <- sum(delhi_diversity$Microcystis == "Present")
n_absent  <- sum(delhi_diversity$Microcystis == "Absent")

# Calculate rank-biserial correlation as the effect size.

rank_sum_MA <- sum(
  rank(delhi_diversity$Shannon)[
    delhi_diversity$Microcystis == "Present"
  ]
)

U_MA <- rank_sum_MA - (n_present * (n_present + 1)) / 2

r_rb_delhi_shannon <- (2 * U_MA / (n_present * n_absent)) - 1

r_rb_delhi_shannon # -0.8148148

# Summarise Shannon diversity within each Microcystis group.

delhi_diversity %>%
  group_by(Microcystis) %>%
  summarise(
    n = n(),
    Median = median(Shannon),
    Minimum = min(Shannon),
    Maximum = max(Shannon)
  )
#Microcystis     n Median Minimum Maximum
#Absent          3   3.48   2.30     3.75
#Present        18   1.73   0.131    2.81

#SPECIES RICHNESS

# Compare species richness between MA and NMA samples
# using the Wilcoxon rank-sum test.

wilcox.test(
  Richness ~ Microcystis,
  data = delhi_diversity,
  exact = FALSE
)
#p = 0.007717

# Calculate rank-biserial correlation as the effect size.

rank_sum_MA <- sum(
  rank(delhi_diversity$Richness)[
    delhi_diversity$Microcystis == "Present"
  ]
)

U_MA <- rank_sum_MA - (n_present * (n_present + 1)) / 2

r_rb_delhi_richness <- (2 * U_MA / (n_present * n_absent)) - 1

r_rb_delhi_richness #-1

# Summarise species richness within each Microcystis group.

delhi_diversity %>%
  group_by(Microcystis) %>%
  summarise(
    n = n(),
    Median = median(Richness),
    Minimum = min(Richness),
    Maximum = max(Richness)
  )

#Microcystis     n Median Minimum Maximum
#Absent          3  321       210     550
#Present        18   72.5      26     187


#MICROBIAL COMMUNITY COMPOSITION COMPARISON IN DELHI

# Extract species abundance data for Delhi samples and
# convert the data into a sample-by-species abundance table.

delhi_abundance <- species %>%
  select(starts_with("D")) %>%
  t() %>%
  as.data.frame()

delhi_abundance$Sample <- rownames(delhi_abundance)

# Add MA/NMA status to the abundance table.

delhi_abundance <- delhi_abundance %>%
  left_join(
    mashanon %>% select(Sample, Microcystis),
    by = "Sample"
  )


#BRAY-CURTIS DISSIMILARITY

# Calculate Bray-Curtis dissimilarity between Delhi samples
# based on their microbial species abundance profiles.

delhi_bray <- vegdist(
  delhi_abundance %>% select(-Sample, -Microcystis),
  method = "bray"
)


#PERMANOVA

# Test whether overall microbial community composition
# differs between MA and NMA samples in Delhi.

delhi_permanova <- adonis2(
  delhi_bray ~ Microcystis,
  data = delhi_abundance,
  permutations = 999
)

delhi_permanova
