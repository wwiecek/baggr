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

  // Selection-adjusted log-likelihood for one observation.
  real sel_loglik_general(real y, real m,
                          real outcome_sd, real z_se,
                          vector c, vector omega) {
    int M = num_elements(c);
    real zabs;
    int j;
    real w_obs;
    real mu_z;
    real sd_z;
    real p;
    real a;

    if (M == 0)
      return normal_lpdf(y | m, outcome_sd);

    zabs = fabs(y) / z_se;

    // Determine observed interval j with [a,b) convention.
    j = 1;
    while (j <= M && zabs >= c[j])
      j += 1;
    w_obs = (j <= M) ? omega[j] : 1.0;

    mu_z = m / z_se;
    sd_z = outcome_sd / z_se;
    p = 0.0;
    a = 0.0;

    for (i in 1:M) {
      real b = c[i];
      p += omega[i] * prob_absz_interval_general(mu_z, sd_z, a, b);
      a = b;
    }
    p += prob_absz_interval_general(mu_z, sd_z,
                                    c[M], positive_infinity());

    return normal_lpdf(y | m, outcome_sd) + log(w_obs) - log(p);
  }

  // Backwards-compatible wrapper for the previous one-scale likelihood.
  real sel_loglik_one(real y, real m, real se,
                      vector c, vector omega) {
    return sel_loglik_general(y, m, se, se, c, omega);
  }
