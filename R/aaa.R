.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is baggr; see vignette('baggr') for tutorial, ?baggr for basic help.")
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM call:")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
}

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}
