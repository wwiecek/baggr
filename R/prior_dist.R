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
#' @param mu    mean of ln(X) for lognormal or location for Student's generalised T
#' @param sigma SD of ln(X) for lognormal or scale for Student's generalised T
#' @param nu    degrees of freedom for Student's generalised T
#' @param shape Shape parameter for LKJ
#' @param order Order of LKJ matrix (typically it does not need to be specified,
#'        as it is inferred directly in the model)
#'
#' @details
#'
#' The prior choice in [baggr::baggr()] is done via distinct arguments for each type of prior,
#' e.g. `prior_hypermean`, or a named list of several passed to `prior`.
#' See the examples below.
#'
#' Notation for priors is "plain-text", in that you can write the distributions as
#' `normal(5,10)`, `uniform(0,100)` etc.
#'
#' Different parameters admit different priors (see [baggr::baggr()] for explanations of
#' what the different `prior_` arguments do):
#'
#' * `prior_hypermean`, `prior_control`, and `prior_beta`
#'    will take `"normal"`, `"uniform"`, `"lognormal"`, and
#'    `"cauchy"` input for scalars.
#'    For a vector hypermean (see `"mutau"` model), it will take any of these
#'    arguments and apply them independently to
#'    each component of the vector, or it can also take a `"multinormal"` argument
#'    (see the example below).
#' * `prior_hypersd`, `prior_control_sd`, and `prior_sigma` will take `"normal"`, `"uniform"`, and `"lognormal"`
#'   but negative parts of the distribution are truncated
#' * `prior_hypercor` allows `"lkj"` input (see Lewandowski _et al._)
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
#' @param to_array for some models where `p` may be larger than 1, Stan will expect
#'                 an array instead of a numeric (even when p == 1), so for compatibiliy
#'                 we return `fam` as an array type
#'
set_prior_val <- function(target, name, prior, p = 1, to_array = FALSE) {
  if(is.null(prior$dist))
    stop("Wrong prior specification")

  if(!(prior$dist %in% names(prior_dist_fam)))
    stop(paste("Prior family must be one of: ",
               paste(names(prior_dist_fam), collapse = ", ")))
  if(prior$dist %in% c("multinormal", "lkj") && p > 1)
    stop("Multi-dimensional priors can't be 'replicated' in set_prior_val")

  if(to_array)
    target[[paste0(name, "_fam")]] <- array(prior_dist_fam[[prior$dist]], dim = 1)
  else
    target[[paste0(name, "_fam")]] <- prior_dist_fam[[prior$dist]]

  if(p > 1)
    target[[paste0(name, "_fam")]] <- rep(target[[paste0(name, "_fam")]], p)

  # For univariates:
  if(!is.null(prior$values)){
    # For now we only allow dimension of 1 (for LKJ) or 3 (for uni, normal, cauchy, t)
    if(length(prior$values) == 2)
      prior$values <- c(prior$values, 0)
    # else if(length(prior$values) > 1)
    # stop("Prior with more than 2 parameters used. Stopping - this is work in progress.")
    if(to_array)
      target[[paste0(name, "_val")]] <- array(prior$values, dim = c(1, length(prior$values)))
    else
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
lognormal <- function(mu, sigma) {
  check_scalar(mu)
  check_scalar(sigma)
  if(sigma <= 0)
    stop("sigma parameter must be positive")
  return(list(dist = "lognormal", values = c(mu, sigma), dimension = 1))
}

#' @rdname priors
#' @export
student_t <- function(nu, mu, sigma) {
  check_scalar(mu)
  check_scalar(sigma)
  if(nu <= 0)
    stop("degrees of freedom (nu) parameter must be positive")
  if(sigma <= 0)
    stop("sigma parameter must be positive")
  return(list(dist = "student_t", values = c(nu, mu, sigma), dimension = 1))
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
