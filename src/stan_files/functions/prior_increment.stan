  real prior_increment (int family, real[] y, real[] pars) {
    real inc;
    if(family == 0)
      inc = uniform_lpdf(y | pars[1], pars[2]);
    if(family == 1)
      inc = normal_lpdf(y | pars[1], pars[2]);
    if(family == 2)
      inc = cauchy_lpdf(y | pars[1], pars[2]);
    if(family == 3)
      inc = student_t_lpdf(y | pars[1], pars[2], pars[3]);
    return inc;
  }

  real prior_increment_vec (int family, real y, vector pars) {
    real inc;
    if(family == 0)
      inc = uniform_lpdf(y | pars[1], pars[2]);
    if(family == 1)
      inc = normal_lpdf(y | pars[1], pars[2]);
    if(family == 2)
      inc = cauchy_lpdf(y | pars[1], pars[2]);
    if(family == 3)
      inc = student_t_lpdf(y | pars[1], pars[2], pars[3]);
    return inc;
  }
