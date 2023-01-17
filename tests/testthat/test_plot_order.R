# unit test for plot.bagger_compare() gg plot ordering

library(baggr)
library(ggplot2)

bg1 <- baggr(schools)
bg2 <- baggr(schools)
bg3 <- baggr(schools)
compare <- baggr_compare(bg1,bg2,bg3)

plot = plot.baggr_compare(compare)

test_that("the big_df order will plot properly",{
	expect_true(is.ggplot(plot))
	expect_true(plot$layers[[1]]$position$reverse)
	expect_true(plot$layers[[2]]$position$reverse)

})
