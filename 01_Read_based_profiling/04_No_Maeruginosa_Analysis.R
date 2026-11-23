
#### Input files

species_file <- "data/species_abundance_matrix.csv"
chemistry_file <- "data/water_chemistry.csv"

species <- read.csv(species_file, check.names = FALSE)
water_chemistry <- read.csv(chemistry_file, check.names = FALSE)


#### Prepare water chemistry data

analysis_df <- water_chemistry %>%
  filter(!is.na(Temperature)) %>%
  select(Sample, Microcystis_abundance, Temperature, Region)


#### Remove M. aeruginosa and rescale abundance

species_no_MA <- species %>%
  filter(Species != "Microcystis aeruginosa") %>%
  mutate(
    across(
      -Species,
      ~ (.x / sum(.x, na.rm = TRUE)) * 100
    )
  )


#### Convert to sample-by-species matrix

community_no_MA <- species_no_MA %>%
  column_to_rownames("Species") %>%
  as.matrix() %>%
  t()


#### Recalculate alpha diversity

diversity_no_MA <- data.frame(
  Sample = rownames(community_no_MA),
  Shannon = diversity(community_no_MA, index = "shannon"),
  Richness = specnumber(community_no_MA)
)


#### NCR-only alpha diversity comparison

ncr_diversity_no_MA <- diversity_no_MA %>%
  filter(grepl("^P", Sample)) %>%
  left_join(
    water_chemistry %>% select(Sample, Microcystis),
    by = "Sample"
  )


#### NCR-only Bray-Curtis PERMANOVA

ncr_community_no_MA <- community_no_MA[
  grepl("^P", rownames(community_no_MA)),
  ,
  drop = FALSE
]

ncr_bray_no_MA <- vegdist(
  ncr_community_no_MA,
  method = "bray"
)

ncr_status_no_MA <- water_chemistry %>%
  filter(grepl("^P", Sample)) %>%
  select(Sample, Microcystis) %>%
  slice(match(rownames(ncr_community_no_MA), Sample))

permanova_ncr_no_MA <- adonis2(
  ncr_bray_no_MA ~ Microcystis,
  data = ncr_status_no_MA,
  permutations = 999
)

permanova_ncr_no_MA


#### Pooled marginal PERMANOVA

community_no_MA_pooled <- community_no_MA[
  analysis_df$Sample,
  ,
  drop = FALSE
]

bray_no_MA_pooled <- vegdist(
  community_no_MA_pooled,
  method = "bray"
)

set.seed(123)

permanova_no_MA_pooled <- adonis2(
  bray_no_MA_pooled ~ Microcystis_abundance + Temperature + Region,
  data = analysis_df,
  permutations = 999,
  by = "margin"
)

permanova_no_MA_pooled
