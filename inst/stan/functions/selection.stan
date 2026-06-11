  // Pr(|Z| in [a,b)) for Z ~ Normal(mu_z, sd_z); b may be +inf.
  real prob_absz_interval_general(real mu_z, real sd_z, real a, real b) {
    real ppos;
    real pneg;

    if (is_inf(b)) {
      ppos = 1.0 - normal_cdf(a, mu_z, sd_z);
      pneg = normal_cdf(-a, mu_z, sd_z);
    } else {
      ppos = normal_cdf(b,  mu_z, sd_z) -
             normal_cdf(a,  mu_z, sd_z);
      pneg = normal_cdf(-a, mu_z, sd_z) -
             normal_cdf(-b, mu_z, sd_z);
    }

    return ppos + pneg;
  }
  
  // Pr( Z in [a,b) ) for Z ~ Normal(mu, 1); b may be +inf, a may be -inf
  real prob_z_interval_general(real mu_z, real sd_z, real a, real b) {
    if (is_inf(b)) {
      return 1.0 - normal_cdf(a, mu_z, sd_z);
    } else if (is_inf(a)) {
      // avoid NaN gradient from normal_cdf(-Inf, ...) in autodiff
      return normal_cdf(b, mu_z, sd_z);
    } else {
      return normal_cdf(b, mu_z, sd_z) - normal_cdf(a, mu_z, sd_z);
    }
  }

  // Selection-adjusted log-likelihood for one observation.
  real sel_loglik_general(real y, real m, real theta,
                          real outcome_sd, real z_se,
                          vector c, vector omega,
                          int symmetric, int possible_selection) {
    int M = num_elements(c);
    real zabs;
    int j;
    real w_obs;
    real mu_z;
    real sd_z;
    real p;
    real a;

    if (M == 0 || possible_selection == 0)
      return normal_lpdf(y | theta, z_se);

    real z = y / z_se;
    if (symmetric)
      z = fabs(z);

    // Determine observed interval j with [a,b) convention.
    j = 1;
    while (j <= M && z >= c[j])
      j += 1;
    w_obs = (j <= M) ? omega[j] : 1.0;

    mu_z = m / z_se;
    sd_z = outcome_sd / z_se;
    p = 0.0;
    a = symmetric ? 0.0 : negative_infinity();

    for (i in 1:M) {
      real b = c[i];
      if (symmetric) {
        p += omega[i] * prob_absz_interval_general(mu_z, sd_z, a, b);
      } else {
        p += omega[i] * prob_z_interval_general(mu_z, sd_z, a, b);
      }
      a = b;
    }
    if (symmetric) {
      p += prob_absz_interval_general(mu_z, sd_z, c[M], positive_infinity());
    } else {
      p += 1.0 - normal_cdf(c[M], mu_z, sd_z);
    }

    return normal_lpdf(y | theta, z_se) + log(w_obs) - log(p);
  }

  // Backwards-compatible wrapper for the previous one-scale likelihood.
  real sel_loglik_one(real y, real m, real se,
                      vector c, vector omega,
                      int symmetric, int possible_selection) {
    return sel_loglik_general(y, m, m, se, se, c, omega,
                              symmetric, possible_selection);
  }
