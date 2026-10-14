#' "Mean and interval" function, including other summaries, calculated for matrix (by column) or vector
#'
#' This function is just a convenient shorthand for getting typical summary statistics.
#'
#' @param y matrix or a vector; for matrices, `mint` is done by-column
#' @param prob Probability mass of the central interval, on the 0--1 scale.
#' A probability of zero returns the median as both endpoints; a probability
#' of one returns the minimum and maximum draws.
#'
#' @param int Compatibility alias for `prob`; supply only one argument.
#' @param digits number of significant digits to [round] values by.
#' @param median return median value?
#' @param sd return SD?
#'
#' @return A named numeric vector for vector input, or a numeric matrix
#'   with one row per input column for matrix input. Statistics are ordered
#'   as `mean`, optional `median` and `sd`, `lower`, `upper`, `width`.
#'   `width` records `prob` even when all draws are unavailable.
#' @export
#'
#' @examples
#' mint(rnorm(100, 12, 5))
#'
mint <- function(y, prob=0.95, digits = NULL, median = FALSE, sd = FALSE,
                 int = NULL){
  prob <- interval_prob(prob, int, missing(prob), missing(int), "int")
  # this is a wrapper that checks for type and applies mintv() accordingly
  if(is.array(y)){
    if(length(dim(y)) == 1)
      y <- as.vector(y)
    if(length(dim(y)) > 2)
      # This behaviour might have to change!
      stop("mint() can't summarise arrays of more than 2 dimensions")
  }
  if (is.matrix(y))
    t(apply(y, 2, mintv, prob = prob, digits = digits, median = median, sd = sd))
  else if (is.numeric(y))
    mintv(y, prob = prob, digits = digits, median = median, sd = sd)
  else
    return(NULL)
}

# Summarise draws while keeping the requested interval width known.
mintv <- function(y, prob = 0.95, digits = NULL, median = FALSE, sd = FALSE) {
  unavailable <- all(is.na(y))
  bounds <- if(unavailable) c(NA_real_, NA_real_) else
    unname(quantile(y, c((1 - prob)/2, (1 + prob)/2)))
  x <- c(mean = if(unavailable) NA_real_ else mean(y))
  if(median)
    x <- c(x, median = if(unavailable) NA_real_ else stats::median(y))
  if(sd)
    x <- c(x, sd = if(unavailable) NA_real_ else
      if(length(y) == 1L) 0 else stats::sd(y))
  x <- c(x, lower = bounds[1], upper = bounds[2])
  if(!is.null(digits))
    x <- round(x, digits = digits)
  c(x, width = prob)
}
