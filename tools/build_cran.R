# Build a CRAN source package without the Windows big-object flag.
# Run from the project root: Rscript tools/build_cran.R
# For submission, use interactive R as described in tools/cran.md.
# Normal development builds retain the flag in configure.win.

library(devtools)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || any(args != "--submit"))
  stop("Usage: Rscript tools/build_cran.R [--submit]")
if (identical(args, "--submit") && !interactive())
  stop("Submission prompts need interactive R; see tools/cran.md.")

local({
  build_dir <- tempfile("baggr-cran-")
  dir.create(build_dir)
  on.exit(unlink(build_dir, recursive = TRUE))

  files <- system2("git", "ls-files", stdout = TRUE)
  stopifnot(length(files) > 0, all(file.exists(files)))
  for (directory in unique(dirname(files)))
    dir.create(file.path(build_dir, directory), recursive = TRUE,
               showWarnings = FALSE)
  stopifnot(all(file.copy(files, file.path(build_dir, files))))

  config <- readLines(file.path(build_dir, "configure.win"))
  stopifnot(sum(grepl("PKG_CXXFLAGS += -Wa,-mbig-obj", config,
                     fixed = TRUE)) == 1L)
  config <- config[!grepl("mbig-obj|Add configuration needed for big object files",
                         config)]
  writeLines(config, file.path(build_dir, "configure.win"))

  if (identical(args, "--submit")) {
    devtools::submit_cran(build_dir)
  } else {
    archive <- devtools::build(build_dir, path = getwd(), manual = FALSE)
    message("CRAN source package: ", archive)
    message("Check this tarball before submitting it; see tools/cran.md.")
  }
})
