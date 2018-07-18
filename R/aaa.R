.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to baggr --development version-- July 2018")
  packageStartupMessage("Use vignette('baggr') for tutorial.")
}

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}
