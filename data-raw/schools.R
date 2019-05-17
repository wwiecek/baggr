schools <- data.frame(
  "group" = paste("School", LETTERS[1:8]),
  "tau" = c(28,8,-3,7,-1,1,18,12),
  "se"  = c(15,10,16,11,9,11,10,18))
devtools::use_data(schools, overwrite = TRUE)
