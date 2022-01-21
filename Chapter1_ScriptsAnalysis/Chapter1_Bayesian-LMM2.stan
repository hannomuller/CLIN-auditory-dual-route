data {
  int<lower=1> N;
  real<lower=0> y[N];
  real x1[N];
  real x2[N];
  real x3[N];
  real x4[N];
  real x5[N];
  int<lower=-1,upper=1> x6[N];
  int<lower=1> J;                            //number of subjects
  int<lower=1> K;                            //number of items
  int<lower=1, upper=J> subj[N];             //subject id
  int<lower=1, upper=K> item[N];             //item id
}

parameters {
  vector[J] u;                                    //subject intercepts
  vector[K] w;                                    //item intercepts
  real<lower=0> sigma_w;                          //item sd
  real<lower=0> sigma_u;                          //subj sd
  real<lower=0> sigma_e;                          //error sd
  real<lower=0> intercept;                        //intercept - must be positive
  real<lower=0, upper=0.99> pp;                    //penalty - must be positive
  real b1;                                        //beta
  real b2;                                       //beta
  real b3;                                       //beta
  real b4;                                        //beta
  real b5;                                       //beta
}

transformed parameters {
  real mu[N];
  real mu_route1[N];
  real mu_route2[N];
  for (i in 1:N) {
    if (x6[i] == -1) {                       // simple cases
      mu[i] = intercept + b1*x1[i] + b2*x2[i] + b3*x3[i] + b5*x5[i] + u[subj[i]] + w[item[i]];
    } else if (x6[i] == 1) {                 // complex cases
      mu_route1[i] = intercept + b1*x1[i] + b2*x2[i] + b3*x3[i] + b5*x5[i] + u[subj[i]] + w[item[i]];
      mu_route2[i] = intercept + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i] + pp + u[subj[i]] + w[item[i]];
      mu[i] = fmin(mu_route1[i], mu_route2[i]);
    }
  }
}

model {
  // priors
  u ~ normal(0, sigma_u);    //subj random effects
  w ~ normal(0, sigma_w);    //item random effects
  pp ~ normal(0.617, 10);   //penalty - high SD due to uncertainty about effect

  // Priors are derived from data in BALDEY that are not in en_words
  intercept ~ normal(7.0098694, 0.21);
  b1 ~ normal(0.4458155, 0.8285539);
  b2 ~ normal(0.2377380, 0.8559615);
  b3 ~ normal(0.0071450, 0.203323);
  b4 ~ normal(-0.0023261, 0.09953084);
  b5 ~ normal(-0.0094699, 0.1296604);

  y ~ normal(mu, sigma_e);
}

generated quantities {
  real y_pred[N];
  for (i in 1:N) {
    y_pred[i] = normal_rng(mu[i], sigma_e);
  }
}