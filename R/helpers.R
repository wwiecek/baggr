paste3_formatter <- function(m, l, u, d){
  paste0(format(m,digits=d),
         " [",
         format(l,digits=d),
         ", ",
         format(u,digits=d),
         "]")
}

extract_logpd_draws <- function(fit) {
  out <- tryCatch(
    rstan::extract(fit, "logpd[1]")[[1]],
    error = function(e) NULL
  )
  if(is.null(out))
    out <- rstan::extract(fit, "logpd")[[1]]
  as.numeric(out)
}
