// generated with brms 2.7.0
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc = K - 1; 
  matrix[N, K - 1] Xc;  // centered version of X 
  vector[K - 1] means_X;  // column means of X before centering 
  for (i in 2:K) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
  } 
} 
parameters { 
  vector[Kc] b;  // population-level effects 
  real temp_Intercept;  // temporary intercept 
} 
transformed parameters { 
} 
model { 
  vector[N] mu = temp_Intercept + Xc * b;
  // priors including all constants 
  target += student_t_lpdf(temp_Intercept | 3, 0, 10); 
  // likelihood including all constants 
  if (!prior_only) { 
    target += bernoulli_logit_lpmf(Y | mu);
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_Intercept = temp_Intercept - dot_product(means_X, b); 
} 
