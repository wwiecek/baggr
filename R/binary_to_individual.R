#' Generate individual-level binary outcome data from an aggregate statistics
#'
#' This is a helper function that is typically used automatically by some of _baggr_ functions,
#' such as when running `model="logit"` in [baggr], when summary-level data are supplied.
#'
#' @param data A data frame with columns `a`, `c` and `b`/`n1`, `d`/`n2`.
#'             (You can also use `ai`, `ci`, `n1i`, `n2i` instead.)
#' @param group Column name storing group
#' @param covariates Column names in `data` that contain group-level variables
#'                   to retain when expanding into individual-level `data.frame`
#' @param rename_group If `TRUE` (default), this will rename the grouping variable
#'                     to `"group"`, making it easier to work with [baggr]
#'
#' See `vignette("baggr_binary")` for an example of use and notation details.
#'
#' @return A data frame with columns `group`, `outcome` and `treatment`.
#' @export
#' @seealso [prepare_ma] uses this function
#'
#' @examples
#' df_yusuf <- read.table(text="
#'   trial        a n1i  c n2i
#'   Balcon      14  56 15  58
#'   Clausen     18  66 19  64
#'   Multicentre 15 100 12  95
#'   Barber      10  52 12  47
#'   Norris      21 226 24 228
#'   Kahler       3  38  6  31
#'   Ledwich      2  20  3  20
#'   ", header=TRUE)
#' bti <- binary_to_individual(df_yusuf, group = "trial")
#' head(bti)
#' # to go back to summary-level data
#' prepare_ma(bti, effect = "logOR")
#' # the last operation is equivalent to simply doing
#' prepare_ma(df_yusuf, group="trial", effect="logOR")
#'
binary_to_individual <- function(data, group = "group",
                                 covariates = c(),
                                 rename_group = TRUE) {
  df_ind <- data.frame()
  if(rename_group)
    group_name <- "group"
  else
    group_name <- group

  if(is.null(data[[group]])){
    message("Automatically created group labels as they were not defined.")
    data[[group]] <- paste("Group", 1:nrow(data))
    # stop("Missing group column")
  }
  if(!is.null(data[["n1i"]]) && is.null(data[["n1"]])) data[["n1"]] <- data[["n1i"]]
  if(!is.null(data[["n2i"]]) && is.null(data[["n2"]])) data[["n2"]] <- data[["n2i"]]
  if(!is.null(data[["ai"]])  && is.null(data[["a"]]))  data[["a"]] <- data[["ai"]]
  if(!is.null(data[["ci"]])  && is.null(data[["c"]]))  data[["c"]] <- data[["ci"]]


  if(is.null(data[["n1"]]) || is.null(data[["n2"]])){
    data <- data[,c(group, "a", "b", "c", "d", covariates)]
    data$n1 <- data$a + data$b
    data$n2 <- data$c + data$d
  }

  data_check <- data.frame(data$a, data$c, data$n1, data$n2)
  non_integer <- lapply(data_check, function(x) any((x %% 1) > 0))
  if(any(unlist(non_integer)))
    stop("Non-integer number of events or non-events present in data")

  for(i in 1:nrow(data)) {
    current_df <- setNames(
      data.frame(
        data[[group]][i],
        c(rep(1, data$n1[i]), rep(0, data$n2[i])),
        c(rep(1, data$a[i]),  rep(0, data$n1[i] - data$a[i]),
          rep(1, data$c[i]),  rep(0, data$n2[i] - data$c[i]))),
      c(group_name, "treatment", "outcome")
    )
    for(cn in covariates)
      current_df[[cn]] <- data[[cn]][i]

    df_ind <- rbind(df_ind, current_df)
  }
  df_ind
}

