## cran-comments for baggr v0.7.11 (May 2025)

Intial four submissions for this version (apologies!!!) had 
(1) some broken/moved permanently URLs and one documentation NOTE. I hope that has now all been fixed, apologies.
(2) NOTE on Windows about -Wa,-mbig-obj flags; I removed them for CRAN release
(3) Accidentally re-added .Rdata in build folder which inflated tarball size

These issues have now been fixed, thank you for your patience.

Test environments, now implemented via GitHub Actions:

* ubuntu-latest, devel, oldrel
* windows-latest
* macos-latest
* generic R-hub configuration (`rhub::rhub_setup()`)

* I also conduct CHECK locally on some Linux distros and Windows 11.

There are no WARNINGs; sometimes the following NOTEs

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
