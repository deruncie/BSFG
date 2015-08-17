// Latent factor model by Rick Farouni
// http://rfarouni.github.io/assets/projects/BayesianFactorAnalysis/BayesianFactorAnalysis.html
// modified by Cody Markelz

data {
  int<lower=1> N;                // number of observations
  int<lower=1> P;                // number of predictors
  matrix[N,P]  Y;                // data matrix of order [N,P]
  int<lower=1> D;                // number of latent dimensions 
}

transformed data {
  int<lower=1> M;
  vector[P] mu;
  M  <- D*(P-D)+ D*(D-1)/2;  // calculate number of non-zero loadings
  mu <- rep_vector(0.0,P);
}

parameters {    
  vector[M] L_t;               // lower triangle elements of L
  vector<lower=0>[D] L_d;      // lower diagonal elements of L
  vector<lower=0>[P] psi;      // vector of variances
  real<lower=0>  mu_psi;
  real<lower=0>  sigma_psi;
  real           mu_lt;
  real<lower=0>  sigma_lt;
}

transformed parameters{
  cholesky_factor_cov[P,D] L;  // lower triangular factor loadings Matrix 
  cov_matrix[P] Q;             // Covariance matrix
{
  int idx1;
  int idx2;
  real zero; 
  zero <- 0;
  for(i in 1:P){
    for(j in (i+1):D){
      idx1 <- idx1 + 1;
      L[i,j] <- zero;       //constrain the upper triangular elements to zero 
    }
  }
  for (j in 1:D) {
      L[j,j] <- L_d[j];     //set diagonal of L
    for (i in (j+1):P) {
      idx2 <- idx2 + 1;
      L[i,j] <- L_t[idx2];  // set values below diagonal of L
    } 
  }
} 
Q <- L*L' + diag_matrix(psi); // recover covariance matrix, add error covariance matrix
}

model {
// the hyperpriors 
   mu_psi ~ cauchy(0, 1);
   sigma_psi ~ cauchy(0,1);
   mu_lt ~ cauchy(0, 1);
   sigma_lt ~ cauchy(0,1);
// the priors 
  L_d ~ cauchy(0,3);
  L_t ~ cauchy(mu_lt,sigma_lt);
  psi ~ cauchy(mu_psi,sigma_psi);
//The likelihood
for( j in 1:N)
    Y[j] ~ multi_normal(mu,Q); 
}


