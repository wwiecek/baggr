# copy files to be compiled into the models folder
copy_models <- function(installation_path = "~/baggr/models_stan/"){
  model_names <- c("joint", "mutau", "rubin")
  for(nm in paste0(model_names, ".stan"))
    file.copy(
      system.file("models", nm, package = "baggr"),
      installation_path
    )
}

.onLoad <- function(libname, pkgname) {
  packageStartupMessage("Welcome to baggr --development version-- May 2018")
  packageStartupMessage("Use vignette('baggr') for tutorial.")

  if(!dir.exists("~/baggr")){
    dir.create("~/baggr")
    dir.create("~/baggr/models_stan")
    packageStartupMessage("Created baggr/ folder in your main directory to store models.
                          Please manually put your models in models_stan/")
    packageStartupMessage("Once it's created try library(baggr) again.")
  }
  if(file.exists("~/baggr/precompiled_models.Rdata")){
    load("~/baggr/precompiled_models.Rdata", envir = .GlobalEnv)
    packageStartupMessage("baggr models loaded from your ~/baggr/ directory.")
    packageStartupMessage("Remember to set options(mc.cores = parallel::detectCores())")

    # logofile <- file("inst/logo.txt")
    # for(i in 1:20)
      # packageStartupMessage(readLines(logofile))
    # close(logofile)

  } else {
    packageStartupMessage("
       First use of package: attempting to compile baggr
       Stan models. This will take a few minutes.
       Expect to see a lot of extra outputs from complier,
       these are not errors.")
    # source('~/baggr/models_stan/model_precompilation.R')
    copy_models()
    packageStartupMessage("Copied .stan models to your baggr/ folder")
    model_precompilation()
    load("~/baggr/precompiled_models.Rdata", envir = .GlobalEnv)
    packageStartupMessage("All baggr Core Stan models compiled in ~/baggr/ folder. Do not remove them.")
  }
}


