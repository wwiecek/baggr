#' Prior distributions in baggr
#'
#' @name priors
#' @description This page provides a list of all available distributions
#' that can be used to specify priors in [baggr()]. These convenience functions
#' are designed only to let user write the priors in the most "natural" way. Apart from
#' passing on the arguments, their only other role is rudimentary check if the distribution
#' is specified correctly.
#'
#' @param location Mean for normal, multivariate normal (in which case `location` is a vector)
#'                 and for Cauchy distributions
#' @param scale SD for Normal, scale for Cauchy
#' @param Sigma Variance-covariance matrix for multivariate normal.
#' @param lower Lower bound for Uniform
#' @param upper Upper bound for Uniform
#' @param shape Shape parameter for LKJ
#' @param order Order of LKJ matrix (typically it does not need to be specified,
#'        as it is inferred directly in the model)
NULL

allowed_prior_dists <- c("uniform", "normal", "multinormal", "cauchy")
check_scalar <- function(x) {
  if(length(x) > 1)
    stop("argument must be scalar")
  if(!is.numeric(x))
    stop("argument must be numeric")
}

set_prior_val <- function(target, name, prior) {
  if(is.null(prior$dist))
    stop("Wrong prior specification")
  if(!(prior$dist %in% allowed_prior_dists))
    stop(paste("Prior family must be one of: ",
               paste(allowed_prior_dists, collapse = ", ")))

  target[[paste0(name, "_fam")]] <- switch(prior$dist,
                                           "uniform" = 0,
                                           "normal" = 1,
                                           "cauchy" = 2,
                                           "multinormal" = 3,
                                           "lkj" = 4)
  # For univariates:
  if(!is.null(prior$values))
    target[[paste0(name, "_val")]] <- prior$values
  # For multivariates:
  if(!is.null(prior$mean))
    target[[paste0(name, "_mean")]] <- prior$mean
  if(!is.null(prior$scale))
    target[[paste0(name, "_scale")]] <- prior$scale

  return(target)
}

#' @rdname priors
#' @export
multinormal <- function(location, Sigma) {
  if(length(location) == 1)
    stop("For 1-dimensional multinormal(), please use normal()")
  x <- try(chol(Sigma), silent = TRUE)
  if(inherits(x, "try-error") || !isSymmetric(Sigma))
    stop("Variance-covariance matrix must be positive semi-definite")
  if(ncol(Sigma) != length(location))
    stop("mean and variance-covariance must have matching dimensions")
  return(list(dist = "multinormal", mean=location, scale=Sigma, dimension = length(location)))
}

#' @rdname priors
#' @export
lkj <- function(shape, order=NULL) {
  check_scalar(shape)
  # order of LKJ is implied by the model, same as in Stan
  return(list(dist = "lkj", values = c(shape), dimension = order))
}

#' @rdname priors
#' @export
normal <- function(location, scale) {
  check_scalar(location)
  check_scalar(scale)
  if(scale <= 0)
    stop("Scale (SD) parameter must be positive")
  return(list(dist = "normal", values = c(location, scale), dimension = 1))
}

#' @rdname priors
#' @export
cauchy <- function(location, scale) {
  check_scalar(location)
  check_scalar(scale)
  if(scale < 0)
    stop("Scale parameter must be positive")
  return(list(dist = "cauchy", values = c(location, scale), dimension = 1))
}

#' @rdname priors
#' @export
uniform <- function(lower, upper) {
  check_scalar(lower)
  check_scalar(upper)
  if(lower > upper)
    stop("In uniform(a,b), a <= b.")
  return(list(dist = "uniform", values = c(lower, upper), dimension = 1))
}
