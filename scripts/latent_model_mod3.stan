// Latent factor model by Rick Farouni
// http://rfarouni.github.io/assets/projects/BayesianFactorAnalysis/BayesianFactorAnalysis.html
// modified by Cody Markelz

data {
  int<lower=1> n;                // number of observations
  int<lower=1> p;                // number of traits
  matrix[p,n]  Y;                // data matrix of order [p,n]
  int<lower=1> K;                // number of latent factors 
}

transformed data {
  int<lower=1> m;
  vector[p] mu;
  vector[K] mu_f;
  m  <- K*(p-K)+ K*(K-1)/2;     // calculate number of non-zero loadings
  mu <- rep_vector(0.0, p);     // mu zero vector
  mu_f <- rep_vector(0.0, K);   // mu_f zero vector
}

parameters {    
  vector[m] L_t;               // lower triangle elements of L
  vector<lower=0>[K] L_f;      // lower diagonal elements of L
  vector<lower=0>[p] psi;      // vector of variances
  real<lower=0>  mu_psi;
  real<lower=0>  sigma_psi;
  real           mu_lt;
  real<lower=0>  sigma_lt;
}

transformed parameters{
  cholesky_factor_cov[p,K] L;            // lower triangular factor loadings Matrix 
  diag_matrix(rep_vector(1, rows(K))) I; // factor identity matrix
  matrix[K,K] L_F;                       // factor matrix 
{
  int idx1;
  int idx2;
  real zero; 
  zero <- 0;
  for(i in 1:p){
    for(j in (i+1):K){
      idx1 <- idx1 + 1;
      L[i,j] <- zero;       // constrain the upper triangular elements to zero 
    }
  }
  for (j in 1:K) {
      L[j,j] <- L_f[j];     // set diagonal of L
    for (i in (j+1):p) {
      idx2 <- idx2 + 1;
      L[i,j] <- L_t[idx2];  // set values below diagonal of L
    } 
  }
} 
// Q <- L*L' + diag_matrix(psi); // recover covariance matrix, add error covariance matrix
// instead calculate the F matrix and put a prior on it
L_mu <- L*F // L [p, f]  
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

for( i in 1:n)
    F[i] ~ multi_normal(mu_f, I); # zero vector

for( i in 1:n)
    Y[i] ~ multi_normal(L_mu, psi); 
}


