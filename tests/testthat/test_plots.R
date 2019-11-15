context("baggr() plots")

library(baggr)
library(ggplot2)

fit <- baggr(schools)
fitplot <- plot(fit)

test_that("Initial themes match expectations",{
  expect_identical(fitplot$theme,
                   baggr_theme_get())
})

# update theme to use a mono font now
baggr_theme_update(text = element_text(family = "mono"))

update_fitplot <- plot(fit)

unaltered_text_elements <- c("face", "colour", "size", "hjust", "vjust", "angle",
                             "lineheight", "margin", "debug")

test_that("Updating plot themes works", {
  expect_identical(update_fitplot$theme$text$family,
               "mono")
  expect_identical(update_fitplot$theme$text[unaltered_text_elements],
                   fitplot$theme$text[unaltered_text_elements])
})

# replace theme info for tex
baggr_theme_replace(text = element_text(family = "mono"))

replace_fitplot <- plot(fit)

# these should all be null since they aren't
# specifically mentioned in the baggr_theme_replace call
testable_elements <- subset(replace_fitplot$theme$text, replace_fitplot$theme$text != "mono" &
                                 !is.logical(replace_fitplot$theme$text))

# this stays false, not a relevant parameter
# force it to null for purposes of testing
testable_elements$inherit.blank <- NULL

test_that("Replacing plot themes works", {
  expect_identical(replace_fitplot$theme$text$family,
                   "mono")
  lapply(testable_elements, expect_null)
})
