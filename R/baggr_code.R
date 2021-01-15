#' Show Stan code for baggr models or objects
#'
#' @param model either a `baggr` object (fitted model) or one of
#'        `"rubin"`, `"mutau"`, `"rubin_full"`
#' @return Nothing is returned in R. Stan code will be opened externally
#'         (e.g. via notepad).
#' @export

show_model <- function(model) {
  if(class(model) == "baggr")
    model_name <- paste0(model$model, ".stan")
  else if(model %in% names(model_names)) #model_names is helper object
    model_name <- paste0(model, ".stan")
  else
    stop("Invalid 'model' argument.")
  model_file <- system.file("models", model_name, package = "baggr")
  file.show(model_file)
}
