# Prepare Imdad et al. vitamin A mortality meta-analysis data.
# The package object is analysis-ready summary data for baggr().

vitamin_a <- read.csv(
  "data-raw/vitamin_a_imdad2022.csv",
  header = FALSE,
  col.names = c("group", "tau", "se"),
  strip.white = TRUE
)

vitamin_a <- subset(vitamin_a, group != "Lin 2008")
row.names(vitamin_a) <- NULL

stopifnot(
  nrow(vitamin_a) == 18,
  all(c("group", "tau", "se") %in% names(vitamin_a)),
  all(is.finite(vitamin_a$tau)),
  all(is.finite(vitamin_a$se)),
  all(vitamin_a$se > 0)
)

save(vitamin_a, file = "data/vitamin_a.rda", compress = "bzip2")
