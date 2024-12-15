
data {
  int<lower=0> N; // Data points
  int<lower=0> S; // Number of participants
  vector[N] y; // acceptability of Adjective
  vector[N] degree; // degree
  array[N] int<lower=1, upper=S> subject; // subject identifier
}

parameters {
  real m_mu;
  real<lower=0> s_mu;
  real m_sig;
  real<lower=0> s_sig;
  real<lower=0> epsilon;
  //vector[S] mu;
  //vector<lower=0>[S] sigma;
  matrix[2,S] z_u;  // normed random effects by participant for mu and sigma
  cholesky_factor_corr[2] L_u; // subj correlation matrix


}

transformed parameters {
  matrix[2,S] u;  //random effects
  vector[S] mu;
  vector<lower=0>[S] sigma;
  vector<lower=0,upper=1>[N] pred;
  u = diag_pre_multiply([s_mu,s_sig]', L_u) * z_u;
  for (s in 1:S) {
    mu[s] = m_mu + u[1,s];
    sigma[s] = exp(m_sig + u[2,s]);
  }
  for (i in 1:N){
    pred[i] = normal_cdf(degree[i]|mu[subject[i]],sigma[subject[i]]);
  }
}

model {
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  m_mu ~ normal(0,1);
  m_sig ~ normal(-2,2);
  s_mu ~ gamma(1.2,1.5);
  s_sig ~ gamma(1.5,4);
  epsilon ~  gamma(1.5,7); 
  y ~ normal(pred, epsilon);
}

