#' Extract prior settings from a baggr model
#'
#' @param bg A fitted [baggr()] model object.
#'
#' @return A named list of priors with class `baggr_prior` that can be
#'   supplied to the `prior` argument of a subsequent [baggr()] call.
#' @export
get_prior <- function(bg){
  if(!inherits(bg, "baggr"))
    stop("`get_prior` expects an object of class 'baggr'.")

  fp <- bg$formatted_prior
  if(is.null(fp))
    stop("No formatted priors were found in the supplied object.")

  prior_entries <- grep("^prior_", names(fp), value = TRUE)
  prior_names <- unique(sub("^prior_(.+)_(fam|val|mean|scale)$", "\\1", prior_entries))

  first_row <- function(x){
    if(is.null(x))
      return(NULL)
    if(is.matrix(x) || is.array(x))
      return(as.numeric(x[1,]))
    as.numeric(x)
  }

  translate_prior_name <- function(name){
    switch(name,
           "sel" = "selection",
           name)
  }

  priors <- lapply(prior_names, function(nm){
    fam_code <- fp[[paste0("prior_", nm, "_fam")]]
    fam_code <- as.numeric(fam_code)[1]
    dist_name <- names(prior_dist_fam)[prior_dist_fam == fam_code]

    if(length(dist_name) == 0)
      stop("Unrecognised prior family in formatted_prior.")

    vals  <- first_row(fp[[paste0("prior_", nm, "_val")]])
    means <- fp[[paste0("prior_", nm, "_mean")]]
    scales <- fp[[paste0("prior_", nm, "_scale")]]

    if(is.matrix(means) || is.array(means))
      means <- as.numeric(means[1,])
    if(is.array(scales) && length(dim(scales)) > 2)
      scales <- scales[, , 1]

    if(length(vals) == 3 && vals[3] == 0)
      vals <- vals[1:2]

    switch(dist_name,
           "normal"     = normal(vals[1], vals[2]),
           "uniform"    = uniform(vals[1], vals[2]),
           "cauchy"     = cauchy(vals[1], vals[2]),
           "lognormal"  = lognormal(vals[1], vals[2]),
           "student_t"  = student_t(vals[1], vals[2], vals[3]),
           "lkj"        = lkj(vals[1]),
           "multinormal"= multinormal(means, scales),
           stop("Unrecognised prior family in formatted_prior."))
  })

  names(priors) <- vapply(prior_names, translate_prior_name, character(1))

  class(priors) <- c("baggr_prior", class(priors))
  priors
}

#' @export
print.baggr_prior <- function(x, ...){
  bullet_points <- paste0("* ", names(x), ": ",
                          vapply(x, print_dist, character(1)))
  cat(paste(bullet_points, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
