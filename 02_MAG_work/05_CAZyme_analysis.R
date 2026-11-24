
library(tidyverse)
library(vegan)

#### Input

all_cazy_tax <- read_csv("data/all_cazy_tax.csv")

#### CAZyme matrix preparation

cazy_matrix <- all_cazy_tax %>%
  select(-MAG, -Sample.ID, -Phylum, -Class, -Order,
         -Family, -Genus, -Species)

taxonomy_filt <- all_cazy_tax

#### Jaccard distance and PCoA

jaccard_dist <- vegdist(cazy_matrix, method = "jaccard", binary = TRUE)

pcoa_res <- cmdscale(jaccard_dist, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Phylum <- taxonomy_filt$Phylum

eig_vals <- pcoa_res$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0]) * 100

pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)

#### PERMANOVA

set.seed(123)

permanova_res <- adonis2(
  jaccard_dist ~ Phylum,
  data = taxonomy_filt,
  permutations = 999,
  by = "margin"
)

permanova_res

R2 <- round(permanova_res$R2[1], 3)
pval <- permanova_res$`Pr(>F)`[1]

p_label <- ifelse(
  pval < 0.001,
  "p < 0.001",
  paste0("p = ", round(pval, 3))
)

perm_text <- paste0("PERMANOVA\nR² = ", R2, "\n", p_label)

#### CAZyme family richness

cazy_matrix_rich <- as.matrix(cazy_matrix)
rownames(cazy_matrix_rich) <- taxonomy_filt$MAG

cazy_richness <- rowSums(cazy_matrix_rich)

richness_df <- data.frame(
  MAG = rownames(cazy_matrix_rich),
  Richness = cazy_richness
) %>%
  left_join(taxonomy_filt, by = "MAG")

phylum_counts <- richness_df %>%
  count(Phylum)

richness_df <- richness_df %>%
  left_join(phylum_counts, by = "Phylum") %>%
  mutate(Phylum_label = paste0(Phylum, "\n(n=", n, ")"))

#### PERMDISP among phyla

phylum_sizes <- taxonomy_filt %>%
  count(Phylum, name = "n_MAGs")

valid_phyla <- phylum_sizes %>%
  filter(n_MAGs >= 5) %>%
  pull(Phylum)

length(valid_phyla)
unique(valid_phyla)

keep_phylum_idx <- taxonomy_filt$Phylum %in% valid_phyla

taxonomy_disp <- taxonomy_filt[keep_phylum_idx, ]

cazy_matrix_disp <- as.matrix(cazy_matrix[keep_phylum_idx, ])
rownames(cazy_matrix_disp) <- taxonomy_disp$MAG

jaccard_disp <- vegdist(
  cazy_matrix_disp,
  method = "jaccard",
  binary = TRUE
)

dispersion <- betadisper(
  jaccard_disp,
  group = taxonomy_disp$Phylum
)

set.seed(123)

dispersion_test <- permutest(
  dispersion,
  permutations = 999
)

dispersion_test

disp_df <- data.frame(
  MAG = names(dispersion$distances),
  Distance_to_centroid = dispersion$distances
) %>%
  left_join(taxonomy_disp, by = "MAG")

summary_df <- disp_df %>%
  group_by(Phylum) %>%
  summarise(
    IQR = IQR(Distance_to_centroid),
    Median = median(Distance_to_centroid),
    SD = sd(Distance_to_centroid),
    .groups = "drop"
  )

#### GH and PL repertoire summaries

gh_cols <- grep("^GH", colnames(all_cazy_tax), value = TRUE)
pl_cols <- grep("^PL", colnames(all_cazy_tax), value = TRUE)

cazy_filtered <- all_cazy_tax %>%
  filter(rowSums(across(all_of(c(gh_cols, pl_cols)))) > 0)

cazy_counts <- cazy_filtered %>%
  mutate(
    GH_sum = rowSums(across(all_of(gh_cols))),
    PL_sum = rowSums(across(all_of(pl_cols)))
  )

write.csv(
  cazy_counts,
  "results/cazy_all_MAGs.csv",
  row.names = FALSE
)

gh_pl_summary <- cazy_counts %>%
  summarise(
    GH_only = sum(GH_sum > 0 & PL_sum == 0),
    GH_PL = sum(GH_sum > 0 & PL_sum > 0)
  )

gh_pl_summary

#### CAZyme family richness and strategy classification

all_cazy_cols <- grep(
  "^(GH|PL)[0-9]+$",
  colnames(cazy_filtered),
  value = TRUE
)

cazy_counts$Family_Richness <- rowSums(
  cazy_counts[, all_cazy_cols] > 0
)

strategy_analysis <- cazy_counts %>%
  mutate(
    Strategy = case_when(
      PL_sum >= 10 ~ "Specialist (10+ PLs)",
      PL_sum >= 5 ~ "Intermediate (5-9 PLs)",
      PL_sum > 0 ~ "Minor PL User (1-4 PLs)",
      TRUE ~ "GH-Only Scavenger"
    )
  )

table(strategy_analysis$Strategy)

strategy_counts <- strategy_analysis %>%
  group_by(Strategy) %>%
  summarise(Count = n(), .groups = "drop")

strategy_counts$Strategy <- factor(
  strategy_counts$Strategy,
  levels = c(
    "GH-Only Scavenger",
    "Minor PL User (1-4 PLs)",
    "Intermediate (5-9 PLs)",
    "Specialist (10+ PLs)"
  ),
  labels = c(
    "Oligosaccharide Scavengers (PL=0)",
    "Minor Degraders (PL 1–4)",
    "Intermediate Specialists (PL 5–9)",
    "Specialists (PL ≥10)"
  )
)

strategy_counts <- strategy_counts %>%
  mutate(Percent = (Count / sum(Count)) * 100)

#### Fisher enrichment analysis across strategies

cazy_cols <- colnames(strategy_analysis)[
  grepl("^(GH|PL)", colnames(strategy_analysis)) &
    !colnames(strategy_analysis) %in% c("GH_sum", "PL_sum")
]

fisher_results <- map_dfr(cazy_cols, function(fam) {
  
  tbl <- table(strategy_analysis$Strategy, strategy_analysis[[fam]])
  
  if (ncol(tbl) < 2) {
    return(NULL)
  }
  
  test <- fisher.test(
    tbl,
    simulate.p.value = TRUE,
    B = 1e5
  )
  
  data.frame(
    Family = fam,
    p_value = test$p.value
  )
})

fisher_results <- fisher_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

#### Family prevalence across strategies

prop_df <- map_dfr(cazy_cols, function(fam) {
  
  df <- strategy_analysis %>%
    group_by(Strategy) %>%
    summarise(
      proportion = mean(.data[[fam]] == 1),
      .groups = "drop"
    )
  
  df$Family <- fam
  df
})

final_df <- prop_df %>%
  left_join(fisher_results, by = "Family")

sig_df <- final_df %>%
  filter(p_adj < 0.05)

#### Strong ecological differentiation

sig_strong <- sig_df %>%
  group_by(Family) %>%
  mutate(
    max_prop = max(proportion),
    min_prop = min(proportion),
    diff = max_prop - min_prop
  ) %>%
  ungroup() %>%
  filter(diff > 0.3)

length(unique(sig_strong$Family))

enrichment_table <- sig_strong %>%
  group_by(Family) %>%
  arrange(desc(proportion), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

write.csv(
  enrichment_table,
  "results/CAZyme_enrichment_table_family.csv",
  row.names = FALSE
)
