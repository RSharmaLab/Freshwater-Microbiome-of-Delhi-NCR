
#### Fisher's exact test for Microcystis presence by region

ma_table <- matrix(
  c(18, 3,
    19, 21),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Region = c("Delhi", "NCR"),
    Microcystis = c("MA", "NMA")
  )
)

fisher.test(ma_table)
