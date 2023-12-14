paste3_formatter <- function(m, l, u, d){
  paste0(format(m,digits=d),
         " [",
         format(l,digits=d),
         ", ",
         format(u,digits=d),
         "]")
}
