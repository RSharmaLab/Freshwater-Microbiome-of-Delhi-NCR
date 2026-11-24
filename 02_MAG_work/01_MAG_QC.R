#### Load libraries

library(tidyverse)


#### Import data

mgqc <- read_csv("data/MAG_QC_df.csv")


#### Assign regions

mgqc <- mgqc %>%
  mutate(
    Group = case_when(
      str_starts(`Sample ID`, "D") ~ "Delhi",
      str_starts(`Sample ID`, "P") ~ "NCR"
    )
  )


#### N50 comparison

mgqc %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    median_N50 = median(N50, na.rm = TRUE),
    .groups = "drop"
  )

wilcox_N50 <- wilcox.test(
  N50 ~ Group,
  data = mgqc,
  exact = FALSE
)

wilcox_N50


#### Completeness comparison

mgqc %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    median_completeness = median(Completeness, na.rm = TRUE),
    .groups = "drop"
  )

wilcox_completeness <- wilcox.test(
  Completeness ~ Group,
  data = mgqc,
  exact = FALSE
)

wilcox_completeness


#### Rank-biserial correlation

rank_biserial <- function(data, variable) {
  
  NCR <- data[[variable]][data$Group == "NCR"]
  Delhi <- data[[variable]][data$Group == "Delhi"]
  
  NCR <- NCR[!is.na(NCR)]
  Delhi <- Delhi[!is.na(Delhi)]
  
  n_NCR <- length(NCR)
  n_Delhi <- length(Delhi)
  
  combined <- c(NCR, Delhi)
  ranks <- rank(combined, ties.method = "average")
  
  rank_sum_NCR <- sum(ranks[seq_len(n_NCR)])
  
  U_NCR <- rank_sum_NCR - (n_NCR * (n_NCR + 1)) / 2
  
  (2 * U_NCR / (n_NCR * n_Delhi)) - 1
}

r_rb_N50 <- rank_biserial(mgqc, "N50")
r_rb_completeness <- rank_biserial(mgqc, "Completeness")

r_rb_N50
r_rb_completeness
