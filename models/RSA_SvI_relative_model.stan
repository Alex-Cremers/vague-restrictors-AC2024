functions {
  matrix logS1(
    matrix utility,
    real lambda,
    real cost_null,
    real cost_atom,
    real cost_only,
    real cost_every,
    real cost_no,
    real cost_adj,
    int quant,
    vector loc_n_list
  ) {
    matrix[256,6] loc_logS1 = rep_matrix(0,256,6);
    loc_logS1[,1] = lambda * (utility[,1] - cost_null);
    loc_logS1[,2] = lambda * (utility[,2] - loc_n_list*cost_atom);
    loc_logS1[,3] = - lambda * (loc_n_list*cost_atom + cost_only); // the log-probability is always 0 with the strong exhasutive answer
    // doing everything conditional on value being finite to avoid having -Inf pseudo-dependent on parameters that mess with gradient.
    real cost_quant = (quant ? cost_every : cost_no);
    for (i in 1:256) {
      loc_logS1[i,4] = is_inf(utility[i,3]) ? negative_infinity() : lambda * (utility[i,3] - cost_quant);
      loc_logS1[i,5] = is_inf(utility[i,4]) ? negative_infinity() : lambda * (utility[i,4] - cost_quant - cost_adj);
      loc_logS1[i,6] = is_inf(utility[i,5]) ? negative_infinity() : lambda * (utility[i,5] - cost_quant - cost_adj);
    }
    for (i in 1:256) {
      loc_logS1[i,] = loc_logS1[i,] - log_sum_exp(loc_logS1[i,]);
    }
    return loc_logS1;
  }
  matrix L1(
    matrix loc_logS1,
    vector loc_log_prior
  ) {
    matrix[256,2] loc_L1;
    loc_L1 = loc_logS1[,5:6];
    for (i in 1:2) {
      loc_L1[,i] = loc_L1[,i] + loc_log_prior;
    }
    loc_L1 = loc_L1 - log_sum_exp(loc_L1);
    return exp(loc_L1);
  }
}

data {
  int<lower=0> N; // number of data points
  int<lower=0> N_measured; // number of data points for which we have a measured posterior
  int N_Id;  // number of HITs
  int N_adj;
  array[N_measured] int measured_indices; // indices of data points for which we have a measured posterior
  array[N_Id] int<lower=1,upper=N_adj> adjective; // adjective for each HIT
  array[N_Id] int<lower=0,upper=1> quantifier; // quantifier for each HIT
  array[N_Id] matrix[256,5] all_utilities; // array of precomputed utility matrices for each message in each world, by participant.
  //(ideally, rewrite this model with a double array of vectors instead of simple array of matrices)
  vector[256] n_list; // number of positive answers in each world (for cost of list messages). Technically an integer but enters multiplications with real
  array[8,128] int<lower=0,upper=256> worlds_indices;
  // Important: the data is assumed to be sorted by Id (so we can easily deal build the prediction vector)
  array[N] int<lower=1,upper=N_Id> Id; // HIT identifier of each data point
  array[N] int<lower=1,upper=8> probe; // position on the degree scale
  vector<lower=0,upper=1>[N_measured] y;  // measured posterior
}


transformed data {
  array[N_Id] row_vector[N_measured] id_index;    // used in the computation of log-lik for LOO
  for (i in 1:N_Id){
    for (n in 1:N_measured){
      id_index[i][n] = (Id[measured_indices[n]]==i);
    }
  }
}

parameters {
  real<lower=0> lambda;
  real<lower=0> cost_null;
  real<lower=0> cost_atom; 
  real<lower=0> cost_only;
  real<lower=0> cost_every;
  real<lower=0> cost_no;
  vector<lower=0>[N_adj] cost_adj;
  real<lower=0> eps;
}

transformed parameters {
  vector[N] prediction;
  vector<lower=0,upper=1>[N_Id] pExh;
  { // use a block to prevent the S1/L1 variables from being saved and creating such a mess.
  array[N_Id] matrix[256,6] all_logS1;
  array[N_Id] matrix[256,2] all_L1;
  for (i in 1:N_Id) {
    all_logS1[i] = logS1(all_utilities[i],lambda,cost_null,cost_atom,cost_only,cost_every,cost_no,cost_adj[adjective[i]],quantifier[i],n_list);
    all_L1[i] = L1(all_logS1[i],all_utilities[i,,1]); // cheating a bit to retrieve the log prior, which is also the utility of the null message
    pExh[i] = sum(all_L1[i][,2])/sum(all_L1[i][,1]+all_L1[i][,2]);
  }
  for (k in 1:N) {
    prediction[k] = sum(all_L1[Id[k],worlds_indices[probe[k]],]);
    //if(is_nan(prediction[k])){print("Problem with row ",k,". Parameters: ",[
    //  lambda,cost_null,cost_atom,cost_only,cost_every,cost_no,cost_adj]);}
  }
  }
}

model {
  cost_null ~ normal(0,5);
  cost_atom ~ normal(0,0.5);
  cost_only ~ normal(0,1);
  cost_every ~ normal(0,1);
  cost_no ~ normal(0,1.5);
  cost_adj ~ normal(0,1);
  lambda ~ gamma(4,1);
  eps ~ gamma(2,6);
  y ~ normal(prediction[measured_indices], eps);
}


generated quantities { // save log_lik by participant for LOO-CV
  array[N_Id] real log_lik;
  {vector[N_measured] pw_log_lik; // pointwise loglik declared as a local variable to save memory
  for (n in 1:N_measured){
    pw_log_lik[n] = normal_lpdf(y[n]|prediction[measured_indices[n]],eps);
  }
  for (i in 1:N_Id){
    log_lik[i] = id_index[i]*pw_log_lik;
  }}
}


