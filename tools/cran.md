# CRAN submission

Run these commands from the project root. Regenerate changed vignettes from
`.Rmd.orig` sources and commit the generated `.Rmd` files and figures together:

```sh
Rscript vignettes/recompute_vigs.R
Rscript tools/build_cran.R
```

`build_cran.R` copies tracked working files to a temporary directory, removes
`-Wa,-mbig-obj` from its `configure.win`, and builds the source tarball in the
project root. Stage new files before building so that they are included.
The working checkout keeps the flag for ordinary development builds. CRAN
previously objected to the explicit Windows assembler flag. This replaces the
manual configure edit formerly kept on the `cranbuild` branch; submit the tarball
from this script rather than a normal development build.

Check the exact tarball, then submit it only after reviewing the results:

```sh
Rscript -e 'devtools::check_built("baggr_0.9.tar.gz", cran = TRUE, remote = TRUE, manual = TRUE)'
R --vanilla --args --submit
```

In the R session started by the last command, run:

```r
source("tools/build_cran.R")
```

`--submit` uses `devtools::submit_cran()` on the temporary flag-free source
copy and asks for confirmation before uploading. It rebuilds the source, so
commit the checked changes and do not edit them between checking and submission.
Replace the tarball filename when the package version changes. Update
`cran-comments.md` with the release version, actual check environments and
results, and explanations for remaining notes. Confirm the submission email
from CRAN. Keep build archives, saved fits and `.RData` files out of Git.
