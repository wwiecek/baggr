# I follow the approach described at https://ropensci.org/blog/2019/12/08/precompute-vignettes/
# To pre-compute everything needed to render the vignettes, which will save time
# on CRAN servers (checks were above 10 mins)
knitr::knit("vignettes/baggr_binary.Rmd.orig", "vignettes/baggr_binary.Rmd")
# then I manually copy figure/ into vignettes/figure/
