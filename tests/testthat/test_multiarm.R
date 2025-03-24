context("baggr() calls with IPD version of Rubin model")
library(baggr)

skip_on_cran()

set.seed(1990)

N  <- 1000
df <- data.frame(
  treatment = factor(sample(c("A", "B", "C"), N, replace = T), levels = c("A", "B", "C")),
  group = sample(paste("Study", 1:10), N, replace = T)
)
df$cl <- sample(1:10, N, replace = T)
df$outcome_cont <- rnorm(N) + (df$treatment == "A")*0.2 + (df$treatment == "B")*0.4 + (df$treatment == "C")*0.6
df$outcome_bin  <- 1*(df$outcome_cont > 0.2)

bg_n <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "none", iter = 150, refresh=0))
bg_p <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "partial", iter = 150, refresh=0))
bg_f <- expect_warning(baggr(df, outcome = "outcome_cont", pooling = "full", iter = 150, refresh=0))

bg_p <- expect_warning(baggr(df, outcome = "outcome_bin", model = "logit",
                             pooling = "partial", iter = 150, refresh=0, cluster = "cl"))
