These are notes for both human and machine maintainers.

- If making edits that will merge to devel, please make additions to NEWS.md
  under a `# devel (work in progress)` heading at the top, so we don't edit
  NEWS.md headings for previous releases.
- NEWS.md descriptions should be very concise, typically one line per PR, or
  a few lines per major update. See historical examples in NEWS.md before adding.
  Err on the side of being short.
- Development is done via GitHub pull requests into devel, typically one per
  issue.
- GitHub Actions (GHA) runs checks for incoming pull requests. For simple changes
  you may choose to not run all of the time-consuming tests. It may also make sense
  to set `Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")` for a lot of development
  checks, which will make things considerably faster. For targetted changes just 
  run targetted tests before handing over to GHA.
