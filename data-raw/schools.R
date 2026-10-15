schools <- data.frame(
  "group" = paste("School", LETTERS[1:8]),
  "yi" = c(28,8,-3,7,-1,1,18,12),
  "sei"  = c(15,10,16,11,9,11,10,18))
usethis::use_data(schools, overwrite = TRUE)
