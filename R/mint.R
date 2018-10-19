# mint = mean and interval function for matrices and vectors
mint <- function(y, int=0.95, digits = NULL){
  if (class(y) == "matrix")
    t(apply(y, 2, mintv, int = int, digits = digits))
  else if (class(y) == "numeric")
    mintv(y, int = int, digits = digits)
  else
    return(NULL)
}

# vectors only
mintv <- function(y, int=0.95, digits = NULL){
  x <- c(quantile(y, (1-int)/2),
         "mean" = mean(y),
         quantile(y, 1 - (1-int)/2))
  if(!is.null(digits))
    x - round(x, digits=digits)
  x
}
