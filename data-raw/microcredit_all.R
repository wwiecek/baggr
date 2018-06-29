library(dplyr)
library(purrr)
# load("data-raw/microcredit_project_data_with_replication_profit.RData")
load("data-raw/microcredit_project_data.RData")
obj <- ls()

nall_na <- function(x) !all(is.na(x))
get_df <- function(st, n) {
  df <- sapply(obj[grep(st, obj)], function(x) {
    y <- get(x)
    if(length(y) == n)
      return(y)
    return(NA) }) %>%
    data.frame() %>%
    select_if(nall_na)

  #remove things like attanassio_
  colnames(df) <- gsub(paste0(st, "_"), "", colnames(df))
  df
}

study <-  c("angelucci", "attanasio", "augsberg", "banerjee", "crepon", "karlan", "tarozzi")
n_subjects <-  c(21523, 961, 1196, 6863, 5498, 1113, 3113)

microcredit_all <- map2(study, n_subjects, get_df)
vars <- c("consumerdurables", "consumption", "expenditures", "profit", "revenues", "temptation", "treatment")
microcredit <- lapply(microcredit_all, function(x) {
  for(nm in vars)
    if(!(nm %in% names(x)))
      x[[nm]] <- NA
    x[, vars]
}) %>% setNames(study) %>% bind_rows(.id = "study")

devtools::use_data(microcredit, overwrite = TRUE)
