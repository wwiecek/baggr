.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to baggr --development version-- winter 2018")
  packageStartupMessage("Use vignette('baggr') for tutorial.")
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM call:")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
}

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}
