# mint = mean and interval function for matrices and vectors
mint <- function(y, int=0.95, digits = NULL){
  if(class(y) == "array" && length(dim(y)) == 1)
    y <- as.vector(y)
  if (class(y) == "matrix")
    t(apply(y, 2, mintv, int = int, digits = digits))
  else if (class(y) == "numeric")
    mintv(y, int = int, digits = digits)
  else
    return(NULL)
}

# vectors only
mintv <- function(y, int=0.95, digits = NULL){
  if(all(is.na(y))){
    y <- 1:10 #spoof, just to get the labels
    y_spoof <- TRUE
  } else
    y_spoof <- FALSE
  x <- c(quantile(y, (1-int)/2),
         "mean" = mean(y),
         quantile(y, 1 - (1-int)/2))
  if(y_spoof)
    x[1:3] <- NA
  if(!is.null(digits))
    x - round(x, digits=digits)
  x
}
