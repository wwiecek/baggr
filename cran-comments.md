## cran-comments for baggr v0.5.0

Run times have for vignettes and examples have been shortened by a factor of 3-4.
If this is too slow for CRAN servers, we will re-write the vignettes to include
precompiled results.

Test environments:

* Ubuntu 16.04.6 LTS x86_64-pc-linux-gnu (64-bit) R 4.0.0 (2020-04-24) via Travis
* win-builder: devel, old release
* Two Windows 10 PCs locally: R 3.5.1, R 4.0.1
* Ubuntu 20.04 R x86_64 3.6.3 locally
* builder.r-hub.io (platforms = NULL)

Check results are unchanged from v0.4.0:
There are no WARNINGs or ERRORs and 2 NOTEs (on some OS's):

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements. 
