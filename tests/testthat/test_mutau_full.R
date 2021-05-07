df <- data.frame()

K <- 10
N <- c(rep(20, 5), rep(20, 5))
bsl <- c(0,0,0,0,0,1,1,1,1,1)
se <- rep(1, 10)
tau <- c(1,1,1,1,1,0,0,0,0,0)

for(i in 1:K){
  df <- rbind(df,
              data.frame(group = i, outcome = bsl[i] + rnorm(N[i], 0, se[i]), treatment = 0),
              data.frame(group = i, outcome = bsl[i] + tau[i] + rnorm(N[i], 0, se[i]), treatment = 1))
}

prepare_ma(df)

# bg1 <- baggr(prepare_ma(df), model = "mutau")
# bg2 <- baggr(df, model = "mutau_full")
# bg2
