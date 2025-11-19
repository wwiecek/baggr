  // Pr( |Z| in [a,b) ) for Z ~ Normal(mu, 1); b may be +inf
  real prob_absz_interval(real mu, real a, real b) {
    real ppos = normal_cdf(b, mu, 1) - normal_cdf(a, mu, 1);
    real pneg;
    if (is_inf(b)) pneg = normal_cdf(-a, mu, 1);
    else           pneg = normal_cdf(-a, mu, 1) - normal_cdf(-b, mu, 1);
    return ppos + pneg;
  }

  // Selection-adjusted log-likelihood for one observation
  real sel_loglik_one(real y, real m, real se,
  vector c,            // length M, ascending > 0
  vector omega) {      // length M, w_{M+1} = 1 (implicit)
  int M = num_elements(c);
  real zabs = fabs(y) / se;

  // determine observed interval j with [a,b) convention
  int j = 1;
  while (j <= M && zabs >= c[j]) j += 1; // if zabs < c[1] -> j=1; ...; else j=M+1
  real w_obs = (j <= M) ? omega[j] : 1.0;

  // total reporting probability under current mean
  real muz = m / se;
  real p = 0.0;
  real a = 0.0;
  for (i in 1:M) {
    real b = c[i];
    p += omega[i] * prob_absz_interval(muz, a, b);
    a = b;
  }
  p += 1.0 * prob_absz_interval(muz, c[M], positive_infinity());

  return normal_lpdf(y | m, se) + log(w_obs) - log(p);
  }
