# Resolve summary-data conventions before conversion and prior construction.
normalise_summary_columns <- function(data, test_data = NULL, model = NULL,
                                      yi = NULL, sei = NULL,
                                      y0i = NULL, se0i = NULL,
                                      group = "group", group_explicit = FALSE,
                                      silent = FALSE) {
  if(!inherits(data, "data.frame"))
    stop("Can't detect input type because it's not a data.frame")
  columns <- list(yi = "yi", sei = "sei", y0i = "y0i",
                  se0i = "se0i", group = group)
  mappings <- list(yi = yi, sei = sei, y0i = y0i, se0i = se0i)
  supplied <- !vapply(mappings, is.null, logical(1))
  messages <- character()
  legacy_type <- NULL
  if(any(supplied)) {
    for(nm in names(mappings)[supplied]) {
      value <- mappings[[nm]]
      if(!is.character(value) || length(value) != 1 ||
         is.na(value) || !nzchar(value))
        stop("`", nm, "` must be a single column name (a string) or NULL.")
    }
    if(supplied[["yi"]] != supplied[["sei"]])
      stop("Supply both `yi` and `sei` column names.")
    if(supplied[["y0i"]] != supplied[["se0i"]] ||
       (any(supplied[c("y0i", "se0i")]) && !supplied[["yi"]]))
      stop("Supply `y0i` and `se0i` together with `yi` and `sei`.")
  }

  # Explicit individual-level models do not interpret summary columns or slab.
  individual <- !is.null(model) &&
    model %in% c("rubin_full", "mutau_full", "logit", "sslab", "quantiles")
  conventions <- list(
    canonical_wide = c(y0i = "y0i", se0i = "se0i", yi = "yi", sei = "sei"),
    canonical = c(yi = "yi", sei = "sei"),
    alias_wide = c(y0i = "y0", se0i = "se0", yi = "y", sei = "se"),
    alias = c(yi = "y", sei = "se"),
    legacy_wide = c(y0i = "mu", se0i = "se.mu", yi = "tau", sei = "se.tau"),
    legacy = c(yi = "tau", sei = "se"))
  present <- vapply(conventions, function(x) all(x %in% names(data)), logical(1))
  full <- present & lengths(conventions) == 4
  mapping <- NULL
  if(!individual && any(supplied)) {
    mapping <- unlist(mappings[supplied])
    # A custom effect mapping can retain recognised baseline columns.
    if(!supplied[["y0i"]] && any(full)) {
      baseline <- conventions[[which(full)[1]]][c("y0i", "se0i")]
      mapping <- c(baseline, mapping)
    }
  } else if(!individual && any(present)) {
    chosen <- which(present)[1]
    mapping <- conventions[[chosen]]
    # Compare all complete recognised conventions, including held-out data.
    for(candidate in conventions[present]) {
      for(nm in intersect(names(mapping), names(candidate))) {
        for(df in list(data, test_data)) {
          if(!is.null(df) && all(c(mapping[[nm]], candidate[[nm]]) %in% names(df)) &&
             !isTRUE(all.equal(df[[mapping[[nm]]]], df[[candidate[[nm]]]])))
            stop("Conflicting summary columns `", mapping[[nm]], "` and `",
                 candidate[[nm]], "`; supply explicit column mappings.")
        }
      }
    }
    if(grepl("^legacy", names(conventions)[chosen])) {
      legacy_type <- if(length(mapping) == 4) "mutau" else "rubin"
      messages <- if(legacy_type == "mutau")
        paste0("Using `mu`, `se.mu`, ", "`tau`, `se.tau` columns.",
               "Starting in baggr 0.9, the preferred names ",
               "are `y0i`, `se0i`, `yi`, and `sei`.") else
        paste0("Using `tau`and `se` column names. ",
               "Starting in baggr 0.9, the preferred names are `yi` and `sei`.")
    } else if(grepl("^alias", names(conventions)[chosen])) {
      messages <- if(length(mapping) == 4)
        paste0("Detected summary-data columns `y0`, `se0`, `y` and `se`; ",
               "preferred names are `y0i`, `se0i`, `yi` and `sei`.") else
        "Detected summary-data columns `y` and `se`; preferred names are `yi` and `sei`."
    }
  }

  if(!individual && !is.null(mapping) && length(mapping) == 2 && any(full))
    mapping <- c(conventions[[which(full)[1]]][c("y0i", "se0i")], mapping)

  if(!is.null(mapping)) {
    for(label in c("data", "test_data")) {
      df <- if(label == "data") data else test_data
      if(is.null(df)) next
      absent <- setdiff(mapping, names(df))
      if(length(absent))
        stop("Missing mapped summary columns in `", label, "`: ",
             paste(absent, collapse = ", "), ".")
      # Subset first so overlapping custom mappings cannot overwrite a source.
      df[names(mapping)] <- df[unname(mapping)]
      if("n.mu" %in% names(df) && !"n0i" %in% names(df))
        df$n0i <- df[["n.mu"]]
      if("n.tau" %in% names(df) && !"ni" %in% names(df))
        df$ni <- df[["n.tau"]]
      if(label == "data") data <- df else test_data <- df
    }
  }
  summary <- !individual && (!is.null(mapping) ||
    (!is.null(model) && model %in% c("rubin", "mutau")) ||
    check_columns_binary(data, stop = FALSE))
  if(summary && !group_explicit && !group %in% names(data) &&
     "slab" %in% names(data)) {
    data[[group]] <- data[["slab"]]
    if(!is.null(test_data)) {
      if(!"slab" %in% names(test_data))
        stop("Missing study-label column `slab` in `test_data`.")
      test_data[[group]] <- test_data[["slab"]]
    }
    # messages <- c(messages,
      # "Using `slab` as study labels. In baggr, the canonical argument is `group`.")
  }
  if(!silent) for(msg in messages) message(msg)
  list(data = data, test_data = test_data, columns = columns,
       legacy = !is.null(legacy_type), legacy_type = legacy_type,
       messages = messages)
}
