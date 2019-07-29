# When working with individual-level data, check that correct columns are used correctly.

check_columns_numeric <- function(data) {
  numcols <- unlist(lapply(data, is.numeric))
  if(sum(!numcols) > 0)
    stop(paste("Column(s)", paste(names(numcols)[numcols], collapse = ","), "are not numeric"))
}

check_columns <- function(data, outcome, group, treatment)  {

  if(!(is.character(outcome) && is.character(group) && is.character(treatment)))
    stop('Arguments "outcome", "group", "treatment" must be of character type')

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
  if(any(is.na(data[[treatment]])))
    stop("Some of treatment values are NA")
  if(any(is.na(data[[outcome]])))
    stop("Some of outcome values are NA")
  if(any(is.na(data[[group]])))
    stop("Some of group values are NA")

  # Treatment has to be dichotomous
  if(!any((data[[treatment]] = 0) | (data[[treatment]] = 1)))
    stop("Treatment column has to have values 0 or 1")

}
