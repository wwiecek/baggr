## cran-comments for baggr v0.6.21

Previous submission of v0.6.20 failed because check times were over 10 minutes.
I'm hoping it's much less after some updates and will work to shorten this further.

Test environments:

* win-builder: devel, old release
* Windows 10 PC locally: R 4.1.2
* Ubuntu 20.04 R x86_64 4.1.3 locally
* builder.r-hub.io (platforms = NULL)

There are no WARNINGs or ERRORs and 2 NOTEs in some environements:

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements. 
  
CRAN checks previously gave 503 URL error on DOI URL's but I cannot replicate 
this issue.
