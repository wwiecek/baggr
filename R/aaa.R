.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is baggr; see vignette('baggr') for tutorial, ?baggr for basic help.")
  packageStartupMessage("From v0.9 (late 2026) notation for inputs and parameters changed. You can still use old defaults.")
  packageStartupMessage("For multicore CPU use: ")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
}

.onLoad <- function(libname, pkgname) {
  invisible(baggr_theme_set())
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}
