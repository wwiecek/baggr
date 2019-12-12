#' Prior distributions in baggr
#'
#' @name priors
#' @description This page provides a list of all available distributions
#' that can be used to specify priors in [baggr()]. These convenience functions
#' are designed to allow the user to write the priors in the most "natural" way when
#' implementing them in baggr. Apart from
#' passing on the arguments, their only other role is to perform a rudimentary check
#' if the distribution is specified correctly.
#'
#' @param location Mean for normal and multivariate normal (in which case `location` is a vector),
#'                 and median for Cauchy distributions
#' @param scale SD for Normal, scale for Cauchy
#' @param Sigma Variance-covariance matrix for multivariate normal.
#' @param lower Lower bound for Uniform
#' @param upper Upper bound for Uniform
#' @param shape Shape parameter for LKJ
#' @param order Order of LKJ matrix (typically it does not need to be specified,
#'        as it is inferred directly in the model)
#'
#' @details
#'
#' The prior choice in [baggr] is always done via 3 distinct arguments: `prior_hypermean`,
#' `prior_hypersd`, and `prior_hypercor`.
#'
#' These respectively refer to the priors on the average of the effects across
#' the groups (hypermean), the standard deviation of the effects across the groups
#' (hypersd), and the correlation in the distribution of parameters across groups
#' when the model allows multivariate shrinkage (say on control group means and effects).
#'
#' Notation for priors is "plain-text", in that you can write the distributions as
#' `normal(5,10)`, `uniform(0,100)` etc.
#' As with any other argument one has the option to simply input the prior directly,
#' e.g. `prior_hypermean = normal(0,1)` , or by creating a named list of custom priors
#' and then inputting the list to the argument `priors`.
#' See the examples below for more.
#'
#' Different parameters admit different priors:
#'
#' * `prior_hypermean` will take `"normal"`, `"uniform"` and `"cauchy"` input for a scalar mean.
#'    For a vector mean, it will take any of these arguments and apply them independently to
#'    each component of the vector, or it can also take a `"multinormal"` argument
#'    (see the example below).
#' * `prior_hypersd` will take `"normal"` and `"uniform"`
#' * `prior_hypercor` allows `"lkj"` input
#'
#'
#' @author Witold Wiecek, Rachael Meager
#'
#' @references
#' Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe.
#' "Generating Random Correlation Matrices Based on Vines and Extended Onion Method."
#' _Journal of Multivariate Analysis_ 100, no. 9 (October 1, 2009): 1989-2001.
#' https://doi.org/10.1016/j.jmva.2009.04.008.
#'
#' @examples
#' # change the priors for 8 schools:
#' baggr(schools, model = "rubin", pooling = "partial",
#'       prior_hypermean = normal(5,5),
#'       prior_hypersd = normal(0,20))
#'
#' \donttest{
#' # passing priors as a list
#' custom_priors <- list(hypercor = lkj(1), hypersd = normal(0,10),
#'                       hypermean = multinormal(c(0,0),matrix(c(10,3,3,10),2,2)))
#' baggr(microcredit_summary_data, model = "mutau",
#'       pooling = "partial", prior = custom_priors)
#' }

NULL

check_scalar <- function(x) {
  if(length(x) > 1)
    stop("argument must be scalar")
  if(!is.numeric(x))
    stop("argument must be numeric")
}

prior_dist_fam <- c("uniform" = 0,
                    "normal" = 1,
                    "cauchy" = 2,
                    "multinormal" = 3,
                    "lkj" = 4)

set_prior_val <- function(target, name, prior) {
  if(is.null(prior$dist))
    stop("Wrong prior specification")
  if(!(prior$dist %in% names(prior_dist_fam)))
    stop(paste("Prior family must be one of: ",
               paste(names(prior_dist_fam), collapse = ", ")))

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
