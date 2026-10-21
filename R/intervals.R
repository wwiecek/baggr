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

# Summarise samples x rows x effects, reducing only a single effect.
summary_draws <- function(x, prob, effects, rows = NULL, aggregate = FALSE,
                          rename_int = FALSE) {
  if(is.null(dim(x)) || length(dim(x)) == 1L)
    x <- array(x, c(length(x), 1, 1))
  else if(is.matrix(x))
    x <- if(aggregate) array(x, c(nrow(x), 1, ncol(x))) else
      array(x, c(dim(x), 1), dimnames = c(dimnames(x), list(NULL)))
  if(!is.null(rows)) dimnames(x)[[2]] <- rows
  dimnames(x)[[3]] <- effects
  stats <- c("mean", "median", "sd", "lower", "upper", "width")
  out <- array(NA_real_, c(dim(x)[2], length(stats), dim(x)[3]),
               dimnames = list(dimnames(x)[[2]], stats, effects))
  for(i in seq_len(dim(x)[3]))
    for(j in seq_len(dim(x)[2]))
      out[j,,i] <- mint(x[,j,i], prob = prob, median = TRUE, sd = TRUE)[stats]
  if(rename_int) {
    out <- out[, stats != "width",, drop = FALSE]
    dimnames(out)[[2]][match(c("lower", "upper"), stats)] <- interval_names(prob)
  }
  if(length(effects) == 1L)
    out <- effect_slice(out, 1)
  out
}

# Select available statistics in the common order.
interval_stats <- function(x, prob) {
  cols <- if(is.null(dim(x))) names(x) else colnames(x)
  keep <- intersect(c("mean", "median", "sd", "lower", "upper", "width"), cols)
  if(is.null(dim(x))) x[keep] else x[, keep, drop = FALSE]
}

# Arrange printed statistics and label bounds with their percentages.
interval_table <- function(x, prob, include_sd = TRUE) {
  if(is.null(dim(x)))
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  x <- interval_stats(x, prob)
  x <- x[, setdiff(colnames(x), c("width", if(!include_sd) "sd")), drop = FALSE]
  colnames(x)[match(c("lower", "upper"), colnames(x))] <- interval_names(prob)
  attr(x, "effect") <- NULL
  x
}

# Extract an effect as a matrix, including when there is only one row.
effect_slice <- function(x, effect) {
  if(is.matrix(x)) {
    if(!(effect %in% c(1, attr(x, "effect"))))
      stop("Unknown effect")
    return(x)
  }
  if(is.character(effect)) effect <- match(effect, dimnames(x)[[3]])
  out <- matrix(x[,,effect], nrow = dim(x)[1], ncol = dim(x)[2],
                dimnames = dimnames(x)[1:2])
  attr(out, "effect") <- dimnames(x)[[3]][effect]
  out
}
