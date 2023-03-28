# test to check for NA values in categorical and continuous covariables

 library(baggr)

 sch1 <- schools
 sch2 <- schools

 sch1$cov <- c(rnorm(7), NA)
 sch2$cov <- c('a','b','c','d','e','f','g', NA)

 test_that("NA values produce errors with messages",{
 	expect_error(baggr(sch1,cov="cov","\n The covariates for these data include an NA value, which will prevent the model from running. \n Please remove the NA values from the covariates and try again."))
 	expect_error(baggr(sch2,cov="cov","\n The covariates for these data include an NA value, which will prevent the model from running. \n Please remove the NA values from the covariates and try again."))

 })