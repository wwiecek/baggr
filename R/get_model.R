# this is a wrapper for evolving (we hope!) capability
# to grab or create Stan models according to request
# and return a model object ready for sampling() calls
get_model <- function(model) {
  if(!is.null(baggr_models[[model]]))
    return(baggr_models[[model]])

  stop(paste0("No '", model, "' model to return.
              Check that you requested an existing model."))
}
