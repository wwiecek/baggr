## cran-comments for baggr v0.7.11 (May 2025)

(Intial two submissions for this version had some broken/moved permanently 
URLs and one documentation NOTE. I hope that has now all been fixed, apologies.
I attempted to also use `remote` and `incoming` check in `devtools::check()`)

Test environments, now implemented via GitHub Actions:

* ubuntu-latest, devel, oldrel1
* windows-latest
* macos-latest
* generic R-hub configuration (`rhub::rhub_setup()`)

* I also conduct CHECK locally on some Linux distros and Windows 11.

There are no WARNINGs; sometimes the following NOTEs

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

On Windows I added the following compilation flags as instructed by Stan devs
-Wa,-mbig-obj
this was done to avoid compilation errors ("file too big") for users of new 
Stan versions. However, it does seem to sometime trigger a NOTE on Windows (in some 
versions of R) for non-standard flags.
