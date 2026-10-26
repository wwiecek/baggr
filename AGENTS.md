These are notes for both human and machine maintainers.

- If making edits that will merge to devel, please make additions to NEWS.md
  under a `# devel (work in progress)` heading at the top, so we don't edit
  NEWS.md headings for previous releases.
- NEWS.md descriptions should be very concise, typically one line per PR, or
  a few lines per major update. ALWAYS check  historical examples in NEWS.md 
  from prior to v0.9 before adding. Err on the side of being short.
- Do not add NEWS items for simple isolated fixes that are not important to the
  end user. An example: fixing one of the tests.
- Development is done via GitHub pull requests into devel, typically one per
  issue.
- GitHub Actions (GHA) runs checks for incoming pull requests. For simple changes
  you may choose to not run all of the time-consuming tests. It may also make sense
  to set `Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")` for a lot of development
  checks, which will make things considerably faster. For targetted changes just 
  run targetted tests before handing over to GHA.
- `cran-comments.md` should be updated by the maintainer before submitting,
  do not edit it without asking
