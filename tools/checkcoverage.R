library(covr)
cc <- covr::package_coverage(type=c("tests"), quiet=FALSE, pre_clean=FALSE)
print(percent_coverage(cc, by = "line"))
# Uses CODECOV_TOKEN from the environment when required by Codecov.
covr::codecov(coverage = cc)
devtools::load_all(quiet = TRUE) #this should re-compile without covr flags if needed?
