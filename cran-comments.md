## cran-comments for baggr v0.6.18

One of previous submissions (v0.5) was rejected due to checktime being over 10 minutes.
I fixed this starting with v0.6 and monitor checktimes in RHub builder, 
currently checks complete in under 600 seconds.
Run times for vignettes and examples have also been shortened in recent versions.

Test environments:

* win-builder: devel, old release
* Windows 10 PC locally: R 4.1.0
* Ubuntu 20.04 R x86_64 4.0.1 locally
* builder.r-hub.io (platforms = NULL)

There are no WARNINGs or ERRORs and 2 NOTEs (in some environements):

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements. 
