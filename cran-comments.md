## cran-comments for baggr v0.8 (Feb 2026)

Test environments, now implemented via GitHub Actions via new version of rhub

* linux
* macos
* windows

* I also conduct CHECK locally on some Linux distros and Windows 11.
* GHA config also includes ubuntu-latest (oldrel-1), ubuntu-latest (release), ubuntu-latest (devel)

There are no WARNINGs; sometimes the following NOTEs/INFOs, depending on OS:

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Previous submission had one problematic URL in a vignette. Fixed by removing the hyperlink.
