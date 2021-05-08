library(covr)
cc <- covr::package_coverage(type=c("tests"), quiet=FALSE, pre_clean=FALSE)
print(percent_coverage(cc, by = "line"))
covr::codecov(coverage = cc, token = "8b28a681-9708-4e69-ac78-8c8f2651925b")
devtools::load_all(quiet = TRUE) #this should re-compile without covr flags if needed?
