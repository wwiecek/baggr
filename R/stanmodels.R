# file generated via rstantools package

# This file is only intended to be used during the installation process
# nocov start
MODELS_HOME <- "src"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "src", getwd())

stan_files <- dir(file.path(MODELS_HOME, "stan_files"),
                  pattern = "stan$", full.names = TRUE)
stanmodels <- lapply(stan_files, function(f) {
  model_cppname <- sub("\\.stan$", "", basename(f))
  stanfit <- rstan::stanc(f, allow_undefined = TRUE,
                          obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  return(do.call(methods::new, args = c(stanfit[-(1:3)], Class = "stanmodel",
                                        mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
}
)
names(stanmodels) <- sub("\\.stan$", "", basename(stan_files))
rm(MODELS_HOME)
# nocov end
