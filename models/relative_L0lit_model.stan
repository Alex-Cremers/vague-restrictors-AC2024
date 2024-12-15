functions {
  vector mix_L0(
    matrix L0_matrix
  ) {
    vector[256] loc_all_L0;
    for (i in 1:256) {
      loc_all_L0[i] = L0_matrix[i,2];
    }
    return loc_all_L0;
  }
}

data {
  int<lower=0> N; // number of data points
  int<lower=0> N_measured; // number of data points for which we have a measured posterior
  int N_Id;  // number of HITs
  array[N_measured] int measured_indices; // indices of data points for which we have a measured posterior
  array[N_Id] int<lower=1,upper=12> adjective; // adjective for each HIT
  array[N_Id] int<lower=0,upper=1> quantifier; // quantifier for each HIT
  array[N_Id] matrix[256,3] all_L0_posteriors; // array of precomputed meanL0 for lit and exh, in each world, by participant.
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
  real<lower=0> eps;
}

transformed parameters {
  vector[N] prediction;
  { // use a block to prevent the L0 variables from being saved and creating such a mess.
  array[N_Id] vector[256] all_L0;
  for (i in 1:N_Id) {
    all_L0[i] = mix_L0(all_L0_posteriors[i]);
  }
  for (k in 1:N) {
    prediction[k] = sum(all_L0[Id[k],worlds_indices[probe[k]]]);
    //if(is_nan(prediction[k])){print("Problem with row ",k,". Parameters: ",[
    //  lambda,cost_null,cost_atom,cost_only,cost_every,cost_no,cost_adj]);}
  }
  }
}

model {
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


