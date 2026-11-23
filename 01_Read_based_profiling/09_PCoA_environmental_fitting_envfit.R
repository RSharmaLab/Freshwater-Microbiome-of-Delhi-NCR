#### Load libraries

library(tidyverse)
library(readxl)
library(vegan)


#### Import data

combined_species <- read_csv("data/species_abund_combined.csv")
combinedhydro <- read_excel("data/Hydrochemical Analysis_Merged.xlsx")


#### Prepare species abundance matrix

stopifnot("Species" %in% names(combined_species))

spp <- combined_species
rownames(spp) <- spp$Species
spp$Species <- NULL

samples_by_species <- as.data.frame(
  t(as.matrix(spp)),
  stringsAsFactors = FALSE
)

samples_by_species[] <- lapply(
  samples_by_species,
  function(x) as.numeric(as.character(x))
)

samples_by_species[is.na(samples_by_species)] <- 0

samples_by_species <- samples_by_species[
  rownames(samples_by_species) != "DHK",
  ,
  drop = FALSE
]

all_zero_species <- colSums(samples_by_species) == 0

samples_by_species <- samples_by_species[
  ,
  !all_zero_species,
  drop = FALSE
]

rs <- rowSums(samples_by_species)

if (median(rs, na.rm = TRUE) > 1.5) {
  samples_by_species <- samples_by_species / 100
}


#### Prepare environmental data

chemistry_numeric <- combinedhydro %>%
  select(Sample, where(is.numeric)) %>%
  filter(Sample != "DHK") %>%
  filter(if_all(-Sample, ~ !is.na(.x)))

common_samples <- rownames(samples_by_species)[
  rownames(samples_by_species) %in% chemistry_numeric$Sample
]

samples_by_species <- samples_by_species[
  common_samples,
  ,
  drop = FALSE
]

chemistry_numeric <- chemistry_numeric[
  match(common_samples, chemistry_numeric$Sample),
]

chemistry_numeric <- chemistry_numeric %>%
  select(-Sample) %>%
  as.data.frame()

rownames(chemistry_numeric) <- common_samples


#### PCoA and Bray-Curtis distance

bray_dist <- vegdist(samples_by_species, method = "bray")

pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- rownames(samples_by_species)

var_explained <- round(
  100 * pcoa_result$eig / sum(pcoa_result$eig),
  2
)

xlab <- paste0("PCoA1 (", var_explained[1], "%)")
ylab <- paste0("PCoA2 (", var_explained[2], "%)")


#### Environmental fitting

set.seed(123)

efit <- envfit(
  pcoa_result$points,
  chemistry_numeric,
  permutations = 999
)

efit


#### BH-adjusted p-values

efit_pvals <- data.frame(
  Variable = names(efit$vectors$pvals),
  p_raw = efit$vectors$pvals,
  q_BH = p.adjust(efit$vectors$pvals, method = "BH")
)


#### Extract and scale significant environmental vectors

efit_vectors <- as.data.frame(efit$vectors$arrows)
efit_vectors$Variable <- rownames(efit_vectors)

arrow_multiplier <- 0.5

efit_vectors$Dim1 <- efit_vectors$Dim1 * arrow_multiplier
efit_vectors$Dim2 <- efit_vectors$Dim2 * arrow_multiplier

signif_vars <- efit_pvals$Variable[efit_pvals$q_BH < 0.05]

efit_vectors_sig <- efit_vectors %>%
  filter(Variable %in% signif_vars) %>%
  mutate(
    Color = "#2ca02c",
    Pval_Label = paste0(
      "q=",
      signif(
        efit_pvals$q_BH[
          match(Variable, efit_pvals$Variable)
        ],
        2
      )
    )
  )