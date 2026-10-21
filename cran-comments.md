## cran-comments for baggr v0.9 (October 2026)

Test environments, implemented via GitHub Actions:

* ubuntu-latest (release, devel, oldrel-1)
* macos-latest (release)
* windows-latest (release)

* Local CHECK on Ubuntu 24.04, R 4.3.3, RStan 2.32.7 and StanHeaders 2.39.1.

The obsolete URLs printed in the binary vignette have been removed at source,
and all three vignettes have been regenerated. The CRAN source package omits
`-Wa,-mbig-obj` from the Windows configuration, as requested for the previous
submission. Ordinary development builds retain the flag.

The local --as-cran check, including --run-donttest, vignette rebuilds and the
PDF manual, has no ERRORs or WARNINGs and four NOTEs:

N checking installed package size ... NOTE
  installed size is 8.6Mb; libs 6.3Mb.

N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

N checking for future file timestamps ... NOTE
  unable to verify current time
  This is specific to the local check environment.

N checking HTML version of manual ... NOTE
  Initially skipped because tidy and V8 were unavailable. After installing
  them, the same HTML validation and math-rendering checks passed separately
  for all 48 manual pages and 9 equations.

The test suite passed 3,600 assertions; eight tests use skip_on_cran().

A separate URL check reports HTTP 403 responses from the publishers for three
DOI links: 10.1002/14651858.CD008524.pub4, 10.1198/004017005000000517 and
10.1257/aer.20180310. These bibliography identifiers have been retained.
