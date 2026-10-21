# Precompute the .Rmd.orig sources so CRAN builds do not fit Stan models.
# Run from the project root with the current baggr version installed:
# Rscript vignettes/recompute_vigs.R          # all vignettes
# Rscript vignettes/recompute_vigs.R baggr    # just the main vignette
# Commit the generated .Rmd files and figures in vignettes/fig/ together.

vignettes <- commandArgs(trailingOnly = TRUE)
if (!length(vignettes))
  vignettes <- c("baggr", "baggr_binary", "baggr_selection")

local({
  old_dir <- setwd("vignettes")
  on.exit(setwd(old_dir))
  for (vignette in vignettes) {
    knitr::knit(paste0(vignette, ".Rmd.orig"), paste0(vignette, ".Rmd"),
                envir = new.env())
  }
})
