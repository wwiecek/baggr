# fully internal function
# given a data.frame try to infer
# what kind of data we're dealing with

detect_input_type <- function(data) {
  # if(class == "baggr_data")
  if(class(data) != "data.frame")
    stop("Can't detect input type because it's not data.frame")

  if("tau" %in% names(data) && "se" %in% names(data)){
    if("treatment" %in% names(data))
      return("pool_narrow") #will need conversion
    else
      return("pool_noctrl_narrow")
  }

  if(!any(is.na(match(c("tau", "mu", "se.mu", "se.tau"),
                      names(data)))))
    return("pool_wide")

  if("site" %in% names(data) &&
     "treatment" %in% names(data) &&
     nrow(data) > length(unique(data$site)))
    return("individual")

  return("unknown")
}
