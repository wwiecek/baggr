## cran-comments for baggr v0.6.4

Previous submission was rejected due to checktime being over 10 minutes.
We removed some of the compiled code and checked again with winbuilder before submission, 
checks completed in 540 seconds.

Run times have for vignettes and examples have been shortened by a factor of 3-4
compared to the previous submission to CRAN.
If this is too slow for CRAN submission, we will re-write the vignettes to include
pre-compiled results.

Test environments:

* Ubuntu 16.04.6 LTS x86_64-pc-linux-gnu (64-bit) R 4.0.2 (2020-06-22) via Travis
* win-builder: devel, old release
* Windows 10 PC locally: R 4.0.5
* Ubuntu 20.04 R x86_64 4.0.1 locally
* builder.r-hub.io (platforms = NULL)

There are no WARNINGs or ERRORs and 2 NOTEs (in some environements):

N checking installed package size ... NOTE
  sub-directories of 1Mb or more: libs
  
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements. 

Some checks have also flagged "potentially" invalid DOI URLs, but 
I checked manually, multiple times, that the URLs work:

* https://doi.org/10.1136/bmj.d549
* https://doi.org/10.1214/11-STS361


