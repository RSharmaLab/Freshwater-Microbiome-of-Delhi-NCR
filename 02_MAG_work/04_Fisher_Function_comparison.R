
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

#### Fisher enrichment analysis

run_fisher <- function(analysis_data) {
  
  fisher_results <- map_dfr(function_cols, function(function_name) {
    
    contingency_table <- table(
      Function = analysis_data[[function_name]],
      Group = analysis_data$Group
    )
    
    if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
      return(NULL)
    }
    
    fisher_test <- fisher.test(contingency_table)
    
    group_summary <- analysis_data %>%
      group_by(Group) %>%
      summarise(
        total_MAGs = n(),
        function_present = sum(.data[[function_name]] == 1, na.rm = TRUE),
        proportion = function_present / total_MAGs,
        .groups = "drop"
      )
    
    MA_prop <- group_summary %>%
      filter(Group == "MA") %>%
      pull(proportion)
    
    NMA_prop <- group_summary %>%
      filter(Group == "NMA") %>%
      pull(proportion)
    
    enriched_group <- group_summary$Group[
      which.max(group_summary$proportion)
    ]
    
    tibble(
      Function = function_name,
      MA_present = ifelse(length(MA_prop) == 0, NA, MA_prop),
      NMA_present = ifelse(length(NMA_prop) == 0, NA, NMA_prop),
      Odds_Ratio = unname(fisher_test$estimate),
      P_value = fisher_test$p.value,
      Enriched_in = enriched_group
    )
  })
  
  fisher_results <- fisher_results %>%
    mutate(
      FDR = p.adjust(P_value, method = "BH"),
      Significant = ifelse(FDR < 0.05, "Yes", "No"),
      Function = str_remove(Function, "^[A-Z]-S-[0-9]+:")
    ) %>%
    arrange(FDR)
  
  significant_functions <- fisher_results %>%
    filter(FDR < 0.05)
  
  list(
    all = fisher_results,
    significant = significant_functions
  )
}

#### All MAGs

fisher_all <- run_fisher(metabolic)

write.csv(
  fisher_all$all,
  "results/Fisher_function_enrichment_all_MAGs.csv",
  row.names = FALSE
)

write.csv(
  fisher_all$significant,
  "results/Fisher_function_enrichment_significant_all_MAGs.csv",
  row.names = FALSE
)

#### NCR MAGs

fisher_NCR <- run_fisher(metabolic_NCR)

write.csv(
  fisher_NCR$all,
  "results/Fisher_function_enrichment_NCR_MAGs.csv",
  row.names = FALSE
)

write.csv(
  fisher_NCR$significant,
  "results/Fisher_function_enrichment_significant_NCR_MAGs.csv",
  row.names = FALSE
)

#### Near-complete NCR MAGs

fisher_NCR_nearcomplete <- run_fisher(metabolic_NCR_nearcomplete)

write.csv(
  fisher_NCR_nearcomplete$all,
  "results/Fisher_function_enrichment_nearcomplete_NCR_MAGs.csv",
  row.names = FALSE
)

write.csv(
  fisher_NCR_nearcomplete$significant,
  "results/Fisher_function_enrichment_significant_nearcomplete_NCR_MAGs.csv",
  row.names = FALSE
)
