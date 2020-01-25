## ------------------------------------------------------------------------
# the ... parameters in the call to baggr() are
# passed to rstan::sampling
# for documentation see ?rstan::sampling
library(baggr)

fit <- baggr(schools, model = "rubin", 
      pooling = "partial",
      prior_hypermean = normal(0, 10),
      prior_hypersd = cauchy(0, 3),
      control = list(
        adapt_delta = 0.9 # 0.8 by default
      ))

fit

## ------------------------------------------------------------------------
# runs for 10,000 iterations per chain instead of 2,000
fit <- baggr(schools, model = "rubin", pooling = "partial",
      prior_hypermean = normal(0,1), prior_hypersd = cauchy(0,2),
      iter = 10000, control = list(
        adapt_delta = 0.95 # like above, to address divergences
      ))

fit

