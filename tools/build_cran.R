# Run interactively from the project root; stage new files first.
# Update cran-comments.md and regenerate changed vignettes before building.
# Only the temporary CRAN copy omits the Windows big-object flag.
library(devtools)

build_dir <- tempfile("baggr-cran-")
dir.create(build_dir)
files <- system2("git", "ls-files", stdout = TRUE)
for (directory in unique(dirname(files)))
  dir.create(file.path(build_dir, directory), recursive = TRUE,
             showWarnings = FALSE)
stopifnot(all(file.copy(files, file.path(build_dir, files))))
config <- readLines(file.path(build_dir, "configure.win"))
config <- config[!grepl("mbig-obj|Add configuration needed", config)]
writeLines(config, file.path(build_dir, "configure.win"))
archive <- devtools::build(build_dir, path = getwd(), manual = FALSE)

# Run these separately: check, review results, then submit the same copy.
devtools::check_built(archive, cran = TRUE, remote = TRUE, manual = TRUE)

devtools::submit_cran(build_dir)
# After confirming CRAN's email: unlink(build_dir, recursive = TRUE)
