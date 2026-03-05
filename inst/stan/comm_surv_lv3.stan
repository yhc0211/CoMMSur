data {
  int<lower=1> N;                 // number of subjects
  int<lower=1> K;                 // number of mediators
  vector[N]     X;                // exposure
  matrix[N, K]  M;                // mediators
  vector[N]     time;             // observed time
  array[N] int<lower=0, upper=1> event; // 1 = event, 0 = censored
}

parameters {
  // Mediation paths
  vector[K] alpha;
  vector[K] beta;

  real<lower=0> sigma_m;
  real<lower=0> lambda;           // Weibull scale (AFT baseline)
  real<lower=0> rho;              // Weibull shape

  // Mixture component parameters
  real<upper=0> mu_neg;           // mean for negative slab (<= 0)
  real<lower=0> mu_pos;           // mean for positive slab (>= 0)
  vector<lower=0>[2] sigma_slab;  // [sigma_neg, sigma_pos]
  real<lower=0> sigma_spike;      // spike sd around 0

  // Joint mixture weights over 9 cells (neg/null/pos × neg/null/pos)
  simplex[9] Pi;

  // Direct effect parameter
  real gamma;                     // Direct effect of X on outcome T
}

transformed parameters {
  vector[3] mu_comp;
  vector[3] sigma_comp;

  // component indices: 1=neg, 2=spike, 3=pos
  mu_comp[1] = mu_neg;
  mu_comp[2] = 0;
  mu_comp[3] = mu_pos;

  sigma_comp[1] = sigma_slab[1];
  sigma_comp[2] = sigma_spike;
  sigma_comp[3] = sigma_slab[2];
}

model {
  /* ── priors  ─────────────────────── */
  sigma_m     ~ cauchy(0, 2.5);
  lambda      ~ cauchy(0, 2.5);
  rho         ~ lognormal(0, 1);
  gamma       ~ normal(0, 5);

  mu_neg      ~ normal(0, 5);
  mu_pos      ~ normal(0, 5);
  sigma_slab  ~ normal(0, 5);
  sigma_spike ~ normal(0, 0.05);

  Pi          ~ dirichlet(rep_vector(1.0, 9));

  // Component grid: (alpha-comp, beta-comp) pairs
  array[9, 2] int states = {
    {1,1}, {1,2}, {1,3},
    {2,1}, {2,2}, {2,3},
    {3,1}, {3,2}, {3,3}
  };

  // Joint spike-and-slab prior over (alpha_k, beta_k) using Pi[j]
  for (k in 1:K) {
    vector[9] log_w;
    for (j in 1:9) {
      int a = states[j,1];  // component for alpha_k
      int b = states[j,2];  // component for beta_k
      log_w[j] = log(Pi[j])
               + normal_lpdf(alpha[k] | mu_comp[a], sigma_comp[a])
               + normal_lpdf(beta[k]  | mu_comp[b], sigma_comp[b]);
    }
    target += log_sum_exp(log_w);
  }

  // Mediator models: M_k ~ Normal(alpha_k * X, sigma_m)
  for (k in 1:K)
    M[, k] ~ normal(alpha[k] * X, sigma_m);

  // Weibull AFT survival with direct effect gamma
  for (i in 1:N) {
    real eta = dot_product(M[i], beta) + gamma * X[i];
    real lambda_i = lambda * exp(eta / rho);  // AFT-consistent scaling
    if (event[i]) {
      target += weibull_lpdf(time[i]  | rho, lambda_i);
    } else {
      target += weibull_lccdf(time[i] | rho, lambda_i);
    }
  }
}

generated quantities {
  // Posterior probabilities over the 9 joint components for each mediator
  matrix[K, 9] joint_comp_probs;

  // Useful summaries
  vector[K] p_alt;      // P(both non-spike) = sum of corners: (neg,neg),(neg,pos),(pos,neg),(pos,pos)
  vector[K] p_negneg;
  vector[K] p_negpos;
  vector[K] p_posneg;
  vector[K] p_pospos;

  vector[K] theta;      // indirect effect alpha_k * beta_k

  array[9, 2] int states = {
    {1,1}, {1,2}, {1,3},
    {2,1}, {2,2}, {2,3},
    {3,1}, {3,2}, {3,3}
  };

  for (k in 1:K) {
    vector[9] log_p;
    for (j in 1:9) {
      int a = states[j,1];
      int b = states[j,2];
      log_p[j] = log(Pi[j])
               + normal_lpdf(alpha[k] | mu_comp[a], sigma_comp[a])
               + normal_lpdf(beta[k]  | mu_comp[b], sigma_comp[b]);
    }
    {
      vector[9] p = softmax(log_p);
      joint_comp_probs[k] = p';

      //  1:(neg,neg), 3:(neg,pos), 7:(pos,neg), 9:(pos,pos)
      p_negneg[k] = p[1];
      p_negpos[k] = p[3];
      p_posneg[k] = p[7];
      p_pospos[k] = p[9];
      p_alt[k]    = p[1] + p[3] + p[7] + p[9];

      theta[k] = alpha[k] * beta[k];
    }
  }
}
