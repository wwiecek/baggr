# Resolve the old argument name and check the interval probability.
interval_prob <- function(prob, legacy = NULL, prob_missing = TRUE,
                          legacy_missing = TRUE, legacy_name = "interval") {
  if(!legacy_missing) {
    if(!prob_missing)
      stop("Supply only one of prob and ", legacy_name)
    prob <- legacy
  }
  if(!is.numeric(prob) || length(prob) != 1L || is.na(prob) ||
     !is.finite(prob) || prob < 0 || prob > 1)
    stop("prob must be a single finite number between 0 and 1")
  prob
}

# Name the two interval endpoints using R's quantile labels.
interval_names <- function(prob) {
  names(stats::quantile(c(0, 1), c((1 - prob)/2, (1 + prob)/2)))
}

# Select summary statistics by name and label the bounds lci and uci.
interval_stats <- function(x, prob) {
  labels <- interval_names(prob)
  cols <- if(is.null(dim(x))) names(x) else colnames(x)
  lower <- if("lci" %in% cols) "lci" else labels[1]
  upper <- if("uci" %in% cols) "uci" else labels[2]
  keep <- c("mean", intersect(c("sd", "median"), cols), lower, upper)
  if(is.null(dim(x))) {
    x <- x[keep]
    names(x) <- c(head(keep, -2), "lci", "uci")
  } else {
    x <- x[, keep, drop = FALSE]
    colnames(x) <- c(head(keep, -2), "lci", "uci")
  }
  x
}

# Arrange printed statistics and label the bounds with their percentages.
interval_table <- function(x, prob, include_sd = TRUE) {
  if(is.null(dim(x)))
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  x <- interval_stats(x, prob)
  if(!include_sd)
    x <- x[, setdiff(colnames(x), "sd"), drop = FALSE]
  colnames(x)[match(c("lci", "uci"), colnames(x))] <- interval_names(prob)
  x
}

# Extract one effect while retaining row and statistic names.
effect_slice <- function(x, effect) {
  matrix(x[,,effect], nrow = dim(x)[1], dimnames = dimnames(x)[1:2])
}
