context("baggr() calls with full model")
library(baggr)
set.seed(1990)

# Generate 8 schools like IPD data
schools_ipd <- data.frame()
N <- c(rep(10, 4), rep(20, 4))
for(i in 1:8){
  x <- rnorm(N[i])
  x <- (x-mean(x))/sd(x)
  x <- x*schools$se[i]*sqrt(N[i])/1.41 + schools$tau[i]

  y <- rnorm(N[i])
  y <- (y-mean(y))/sd(y)
  y <- y*schools$se[i]*sqrt(N[i])/1.41

  schools_ipd <- rbind(schools_ipd,
                       data.frame(group = schools$group[i], outcome = x, treatment = 1),
                       # This is just so that we don't trip off prepare_ma:
                       data.frame(group = schools$group[i], outcome = y, treatment = 0))
}


bg_n <- expect_warning(baggr(schools_ipd, pooling = "none", iter = 200, refresh=0))
bg_p <- expect_warning(baggr(schools_ipd, pooling = "partial", iter = 200, refresh=0))
bg_f <- expect_warning(baggr(schools_ipd, pooling = "full", iter = 200, refresh=0))


test_that("Different pooling methods work for the full model", {
  expect_is(bg_n, "baggr")
  expect_is(bg_p, "baggr")
  expect_is(bg_f, "baggr")
})

test_that("Basic operations on full data model", {
  expect_error(baggr(schools_ipd, rubbish = 41))
  expect_is(pooling(bg_p)[,,1], "matrix")
  expect_is(plot(bg_p), "gg")
  expect_is(effect_plot(bg_p), "gg")
  expect_is(forest_plot(bg_p), "vpPath")
  bgc <- try(baggr_compare(bg_n, bg_p, bg_f))
  expect_is(bgc, "baggr_compare")
})

test_that("Full model crashes with nonsense inputs", {
  expect_error(baggr(schools_ipd, outcome = 2), "Arguments")
  expect_error(baggr(schools_ipd, group = 2), "Arguments")
  expect_error(baggr(schools_ipd, treatment = 2), "Arguments")
  expect_error(baggr(schools_ipd, treatment = "wrong"), "no column")
  schools_ipd2 <- schools_ipd; schools_ipd2$treatment <- as.character(schools_ipd$treatment)
  expect_error(baggr(schools_ipd2), "has to be numeric")
  schools_ipd2 <- schools_ipd; schools_ipd2$outcome <- as.character(schools_ipd$outcome)
  expect_error(baggr(schools_ipd2), "has to be numeric")

})
