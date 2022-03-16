  real prior_increment_real (int family, real y, vector pars) {
    real inc;
    if(family == 0)
      inc = uniform_lpdf(y | pars[1], pars[2]);
    if(family == 1)
      inc = normal_lpdf(y | pars[1], pars[2]);
    if(family == 2)
      inc = cauchy_lpdf(y | pars[1], pars[2]);
    //3 and 4 are multinormal and lkj, which for now are not incremented
    //by using this function
    if(family == 5)
      inc = lognormal_lpdf(y | pars[1], pars[2]);
    if(family == 6)
      inc = student_t_lpdf(y | pars[1], pars[2], pars[3]);
    return inc;
  }

  real prior_increment_vec (int family, vector y, vector pars) {
    real inc;
    if(family == 0)
      inc = uniform_lpdf(y | pars[1], pars[2]);
    if(family == 1)
      inc = normal_lpdf(y | pars[1], pars[2]);
    if(family == 2)
      inc = cauchy_lpdf(y | pars[1], pars[2]);
    if(family == 5)
      inc = lognormal_lpdf(y | pars[1], pars[2]);
    if(family == 6)
      inc = student_t_lpdf(y | pars[1], pars[2], pars[3]);
    return inc;
  }
