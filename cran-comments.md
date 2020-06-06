## cran-comments for baggr v0.5.0

Test environments:
* Ubuntu 16.04, R 3.6.1 x86_64-pc-linux-gnu via Travis
* win-builder: devel, old release
* Two Windows 10 PCs locally: R 3.5.1, R 3.6.1
* builder.r-hub.io (platforms = NULL)

Check results are unchanged from v0.4.0:
There are no WARNINGs or ERRORs and 2 NOTEs:

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements. 

In R-hub tests 
* Debian Linux, R-devel, GCC ASAN/UBSAN failed due to lack of availability of various dependencies (loo, forestplot, rstan).
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit cannot compile .stan files.
