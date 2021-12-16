# When working with individual-level data, check that correct columns are used correctly.

check_columns_numeric <- function(data) {
  numcols <- unlist(lapply(data, is.numeric))
  if(sum(!numcols) > 0)
    stop(paste("Column(s)", paste(names(numcols)[numcols], collapse = ","), "are not numeric"))
}

is_binary <- function(v) {
  if(identical(as.numeric(sort(unique(v))), c(0,1)))
    return(TRUE)
  else
    return(FALSE)
}

check_columns <- function(data, outcome, group, treatment, stop.for.na = TRUE)  {

  if(!(is.character(outcome) && is.character(group) && is.character(treatment)))
    stop('Arguments "outcome", "group", "treatment" must be of type "character"')

  # Do columns exist?
  if(is.null(data[[outcome]]))
    stop(paste0("There's no column '", outcome, "' in data"))
  if(is.null(data[[group]]))
    stop(paste0("There's no column '", group, "' in data"))
  if(is.null(data[[treatment]]))
    stop(paste0("There's no column '", treatment, "' in data"))

  # Are they the correct type
  if(!(is.numeric(data[[outcome]])))
    stop(paste0("Outcome variable in baggr has to be numeric (continuous)"))
  if(!(is.numeric(data[[treatment]])))
    stop(paste0("Treatment variable in baggr has to be numeric"))

  # NA not allowed
  if(stop.for.na){
    if(any(is.na(data[[treatment]])))
      stop("Some of treatment values are NA")
    if(any(is.na(data[[outcome]])))
      stop("Some of outcome values are NA")
    if(any(is.na(data[[group]])))
      stop("Some of group values are NA")
  }

  # Treatment has to be dichotomous
  if(!any((data[[treatment]] = 0) | (data[[treatment]] = 1)))
    stop("Treatment column has to have values 0 or 1")

}

check_columns_binary <- function(data, stop=TRUE) {
  if(is.null(data$a) || is.null(data$c) ||
     ((is.null(data$b) || is.null(data$d)) && (is.null(data$n1) || is.null(data$n2)))){
    if(stop)
      stop("Binary data must have columns 'a', 'c' and 'b'/'n1', 'd'/'n2'")
    return(0)
  }
  return(1)
}
