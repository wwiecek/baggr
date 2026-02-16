# Agent Instructions

- For future R package compilations via `devtools::load_all()` or `pkgload::load_all()`, always set `quiet = TRUE` unless the user explicitly asks for verbose compiler output.
- Before recompiling or calling `devtools::load_all()` / `pkgload::load_all()`, run `Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")`.
- All fixes must be implemented on new branches created from `devel`.
