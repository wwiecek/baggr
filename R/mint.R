#' "Mean and interval" function, including other summaries, calculated for matrix (by column) or vector
#'
#' This function is just a convenient shorthand for getting typical summary statistics.
#'
#' @param y matrix or a vector; for matrices, `mint` is done by-column
#' @param int probability interval (default is 95 percent) to calculate
#' @param digits number of significant digits to [round] values by.
#' @param median return median value?
#' @param sd return SD?
#'
#' @export
#'
#' @examples
#' mint(rnorm(100, 12, 5))
#'
mint <- function(y, int=0.95, digits = NULL, median = FALSE, sd = FALSE){
  # this is a wrapper that checks for type and applies mintv() accordingly
  if(is.array(y)){
    if(length(dim(y)) == 1)
      y <- as.vector(y)
    if(length(dim(y)) > 2)
      # This behaviour might have to change!
      stop("mint() can't summarise arrays of more than 2 dimensions")
  }
  if (is.matrix(y))
    t(apply(y, 2, mintv, int = int, digits = digits, median = median, sd = sd))
  else if (inherits(y, "numeric"))
    mintv(y, int = int, digits = digits, median = median, sd = sd)
  else
    return(NULL)
}

# mint() applied to vectors
mintv <- function(y, int=0.95, digits = NULL, median = FALSE, sd = FALSE){
  if(all(is.na(y))){
    y <- 1:10 #spoof, just to get the labels
    y_spoof <- TRUE
  } else
    y_spoof <- FALSE
  x <- c(quantile(y, (1-int)/2),
         "mean" = mean(y),
         quantile(y, 1 - (1-int)/2))
  if(median)
    x <- c(x, "median" = median(y))
  if(sd)
    x <- c(x, "sd" = sd(y))
  if(y_spoof)
    x[1:length(x)] <- NA
  if(!is.null(digits))
    x - round(x, digits=digits)
  x
}
