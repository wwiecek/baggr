# fully internal function
# given a data.frame try to infer
# what kind of meta-analysis data we're dealing with

detect_input_type <- function(data,
                              group="group",
                              treatment="treatment",
                              outcome="outcome") {
  # if(class == "baggr_data")
  if(!inherits(data, "data.frame"))
    stop("Can't detect input type because it's not a data.frame")

  # Full baseline/effect sets must precede their effect-only subsets.
  for(cols in list(c("y0i", "se0i", "yi", "sei"),
                   c("y0", "se0", "y", "se"),
                   c("mu", "se.mu", "tau", "se.tau")))
    if(all(cols %in% names(data))) return("pool_wide")
  for(cols in list(c("yi", "sei"), c("y", "se"), c("tau", "se")))
    if(all(cols %in% names(data))) return("pool_noctrl_narrow")
  if(check_columns_binary(data, stop = FALSE))
    return("pool_binary")

  # Individual-level data -----
  if(!is.null(data[[group]])){
    if(nrow(data) > length(unique(data[[group]]))){
      if(is.binary(data[[outcome]]))
        return("individual_binary")
      else
        return("individual")
    }
  }

  # Shrug emoji -----
  return("unknown")
}
