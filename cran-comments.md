## cran-comments for baggr v0.9 (October 2026)

Resubmitting. The first submission of this v0.9 had problems

(1) one of the tests that may have failed stochastically in some v0.8.2 checks
    has been fixed
(2) there were many new gcc-san RcppParallel-related messages in some environments, 
    I do not think they have anything to do with this package and I do not know 
    how I could address that;
    gcc-UBSAN 0.8.2 check dated 2026-10-10 reports Status: OK.

*****

Test environments, implemented via GitHub Actions:

* ubuntu-latest (release, devel, oldrel-1)
* macos-latest (release)
* windows-latest (release)
* for 0.9 added ubuntu-24.04-arm (release) 

* Local CHECK on Ubuntu 24.04, R 4.3.3, RStan 2.32.7 and StanHeaders 2.39.1.

Checks --as-cran pass with two notes:

N checking installed package size ... NOTE
  installed size is 8.6Mb; libs 6.3Mb.

N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
