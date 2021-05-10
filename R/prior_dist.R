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
#'
#' @examples
#' # (these are not the recommended priors -- for syntax illustration only)
#'
#' # change the priors for 8 schools:
#' baggr(schools, model = "rubin", pooling = "partial",
#'       prior_hypermean = normal(5,5),
#'       prior_hypersd = normal(0,20))
#'
#' \donttest{
#' # passing priors as a list
#' custom_priors <- list(hypercor = lkj(1), hypersd = normal(0,10),
#'                       hypermean = multinormal(c(0,0),matrix(c(10,3,3,10),2,2)))
#' microcredit_summary_data <- prepare_ma(microcredit, outcome = "consumption")
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

#' Output a distribution as a string
#'
#' Used for printing nicely formatted outputs when reporting results etc.
#' @param dist distribution name, one of [priors]
#' @return Character string like `normal(0, 10^2)`.
print_dist <- function(dist){
  if(dist$dist == "multinormal")
    return("multinormal(...)")

  paste0(dist$dist, "(",
         paste0(format(dist$values, digits = 2, trim = TRUE), collapse = ", "),
         ifelse(dist$dist == "normal", "^2", ""), ")")
}


#' Add prior values to Stan input for baggr
#'
#' @param target list object (Stan input) to which prior will be added
#' @param name prior name, like `hypermean`, `hypersd`, `hypercor`
#' @param prior one of prior distributions allowed by baggr like [normal]
#' @param p number of repeats of the prior, i.e. when P i.i.d. priors are set for
#'                  P dimensional parameter as in "mu & tau" type of model
#'
set_prior_val <- function(target, name, prior, p = 1) {
  if(is.null(prior$dist))
    stop("Wrong prior specification")

  if(!(prior$dist %in% names(prior_dist_fam)))
    stop(paste("Prior family must be one of: ",
               paste(names(prior_dist_fam), collapse = ", ")))
  if(prior$dist %in% c("multinormal", "lkj") && p > 1)
    stop("Multi-dimensional priors can't be 'replicated' in set_prior_val")

  target[[paste0(name, "_fam")]] <- switch(prior$dist,
                                           "uniform" = 0,
                                           "normal" = 1,
                                           "cauchy" = 2,
                                           "multinormal" = 3,
                                           "lkj" = 4)
  if(p > 1)
    target[[paste0(name, "_fam")]] <- rep(target[[paste0(name, "_fam")]], p)

  # For univariates:
  if(!is.null(prior$values)){
    # For now we only allow dimension of 1 (for LKJ) or 3 (for uni, normal, cauchy, t)
    if(length(prior$values) == 2)
      prior$values <- c(prior$values, 0)
    else if(length(prior$values) > 1)
      stop("Prior with more than 2 parameters used. Stopping - this is work in progress.")
    target[[paste0(name, "_val")]] <- prior$values

    if(p > 1)
      target[[paste0(name, "_val")]] <- t(replicate(p, prior$values))
  }

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
