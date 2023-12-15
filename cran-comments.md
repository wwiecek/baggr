## cran-comments for baggr v0.7.7

Test environments:

* win-builder
* Windows 10 PC locally: R 4.2.2
* Ubuntu 22.04 locally
* builder.r-hub.io (via rhub::check_for_cran)

There are no WARNINGs or ERRORs and 2 NOTEs in some environements:

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Some Rhub environments gave preperrors, 
but that does not seem related to the package.

On Windows I added the following compilation flags as instructed by Stan devs
-Wa,-mbig-obj
this was done to avoid compilation errors ("file too big") for users of new 
Stan versions. However, it does seem to trigger a NOTE on Windows (in some 
versions of R at least)

