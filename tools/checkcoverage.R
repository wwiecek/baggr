library(covr)
cc <- covr::package_coverage(type=c("tests"), quiet=FALSE, pre_clean=FALSE)
print(percent_coverage(cc, by = "line"))
covr::codecov(coverage = cc)
