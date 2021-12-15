#' Generate individual-level binary outcome data from an aggregate statistics
#'
#' This is a helper function that is typically used automatically by some of _baggr_ functions,
#' such as when running `model="logit"` in [baggr], when summary-level data are supplied.
#'
#' @param data A data frame with columns `a`, `c` and `b`/`n1`, `d`/`n2`.
#'
#' See `vignette("baggr_binary")` for an example of use and notation details.
#' @param group Column name storing group
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
binary_to_individual <- function(data, group = "group") {
  df_ind <- data.frame()
  if(is.null(data[[group]]))
    stop("Missing group column")

  if(is.null(data$n1) || is.null(data$n2)){
    data <- data[,c(group, "a", "b", "c", "d")]
    data$n1 <- data$a + data$b
    data$n2 <- data$c + data$d
  }

  data_check <- data.frame(data$a, data$c, data$n1, data$n2)
  non_integer <- lapply(data_check, function(x) any((x %% 1) > 0))
  if(any(unlist(non_integer)))
    stop("Non-integer number of events or non-events present in data")

  for(i in 1:nrow(data)) {
    df_ind <-
      rbind(df_ind,
            data.frame(
              group     = data[[group]][i],
              treatment = c(rep(1, data$n1[i]), rep(0, data$n2[i])),
              outcome   = c(rep(1, data$a[i]),  rep(0, data$n1[i] - data$a[i]),
                            rep(1, data$c[i]),  rep(0, data$n2[i] - data$c[i]))))
  }
  df_ind
}

