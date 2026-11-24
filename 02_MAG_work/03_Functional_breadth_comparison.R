
library(tidyverse)

#### Input and preparation

metabolic <- read_csv("data/biogeochemical_cycling_df.csv")
mag_qc <- read_csv("MAG_QC/data/MAG_QC_df.csv")

metabolic <- metabolic %>%
  left_join(
    mag_qc %>% select(MAG_ID, Completeness, Contamination),
    by = "MAG_ID"
  )

mawaterbodies <- c(
  "DRH", "DRG", "DAP", "DBH", "DDW", "DHK", "DHS", "DJW",
  "DJP", "DLG", "DNZ", "DNL", "DPN", "DPQ", "DRS", "DSD",
  "DSL", "DSP", "PAJ", "PBS", "PBV", "PDA", "PDE", "PDH",
  "PDN", "PDU", "PGA", "PHK", "PHL", "PIK", "PMA", "PMP",
  "PNN", "PPA", "PPP", "PRO", "PSE"
)

metabolic <- metabolic %>%
  mutate(Group = ifelse(`Sample ID` %in% mawaterbodies, "MA", "NMA"))

function_cols <- colnames(metabolic)[9:40]

metabolic_NCR <- metabolic %>%
  filter(str_starts(`Sample ID`, "P"))

metabolic_NCR_nearcomplete <- metabolic_NCR %>%
  filter(Completeness > 90, Contamination < 5)

#### Functional breadth comparison

run_breadth <- function(analysis_data) {
  
  functional_breadth_df <- analysis_data %>%
    mutate(
      Functional_Breadth = rowSums(
        select(., all_of(function_cols)) == 1,
        na.rm = TRUE
      )
    )
  
  breadth_summary <- functional_breadth_df %>%
    group_by(Group) %>%
    summarise(
      MAGs = n(),
      Median = median(Functional_Breadth, na.rm = TRUE),
      Minimum = min(Functional_Breadth, na.rm = TRUE),
      Maximum = max(Functional_Breadth, na.rm = TRUE),
      .groups = "drop"
    )
  
  wilcox_result <- wilcox.test(
    Functional_Breadth ~ Group,
    data = functional_breadth_df,
    exact = FALSE
  )
  
  n_MA <- sum(functional_breadth_df$Group == "MA")
  n_NMA <- sum(functional_breadth_df$Group == "NMA")
  
  rank_sum_MA <- sum(
    rank(functional_breadth_df$Functional_Breadth)[
      functional_breadth_df$Group == "MA"
    ]
  )
  
  U_MA <- rank_sum_MA - (n_MA * (n_MA + 1)) / 2
  
  r_rb <- (2 * U_MA / (n_MA * n_NMA)) - 1
  
  print(breadth_summary)
  print(wilcox_result)
  
  cat("MA MAGs:", n_MA, "\n")
  cat("NMA MAGs:", n_NMA, "\n")
  cat("Wilcoxon W:", as.numeric(wilcox_result$statistic), "\n")
  cat("Wilcoxon p-value:", wilcox_result$p.value, "\n")
  cat("Rank-biserial correlation:", r_rb, "\n")
  
  invisible(list(
    data = functional_breadth_df,
    summary = breadth_summary,
    wilcox = wilcox_result,
    rank_biserial = r_rb
  ))
}

#### All MAGs

breadth_all <- run_breadth(metabolic)

#### NCR MAGs

breadth_NCR <- run_breadth(metabolic_NCR)

#### Near-complete NCR MAGs

breadth_NCR_nearcomplete <- run_breadth(metabolic_NCR_nearcomplete)
