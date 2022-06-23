  real realprior_lpdf (real theta, int family, vector pars) {
    if(family == 0) return uniform_lpdf(theta | pars[1], pars[2]);
    else if(family == 1) return normal_lpdf(theta | pars[1], pars[2]);
    else if(family == 2) return cauchy_lpdf(theta | pars[1], pars[2]);
    //3 and 4 are multinormal and lkj, which for now are not incremented
    //by using this function
    else if(family == 5) return lognormal_lpdf(theta | pars[1], pars[2]);
    // if(family == 6)
    else return student_t_lpdf(theta | pars[1], pars[2], pars[3]);
  }

  real vecprior_lpdf (vector theta, int family, vector pars) {
    if(family == 0) return uniform_lpdf(theta | pars[1], pars[2]);
    else if(family == 1) return normal_lpdf(theta | pars[1], pars[2]);
    else if(family == 2) return cauchy_lpdf(theta | pars[1], pars[2]);
    else if(family == 5) return lognormal_lpdf(theta | pars[1], pars[2]);
    // if(family == 6)
    else return student_t_lpdf(theta | pars[1], pars[2], pars[3]);
  }
