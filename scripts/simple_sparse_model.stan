// simple sparse factor model
// latent portion
// Cody Markelz -- markelz@gmail.com

data {
  int<lower=1> n;                // number of individuals 1 ...n
  int<lower=1> p;                // number of traits 1 ... p
  int<lower=1> K;                // number of latent factors 
  matrix[n,p]  Y;                // data matrix of order [n,p]
  matrix[K,K]  I;
}

transformed data {
  vector[K] mu_f;
  mu_f <- rep_vector(0.0, K);   // mu_f zero vector
}

parameters{
  vector[K] F[n];
  real<lower=0>   sigma;
  real            psi;
  matrix[p,K] Lambda;
  real            lam;
}

transformed parameters{
//  diag_matrix(rep_vector(1, K)) I; // factor identity matrix
}

model {
  psi ~ gamma(3/2, 3/2);
  sigma ~ cauchy(0, 1);
  
  for( j in 1:p) {
    for (k in 1:K) {
      Lambda[j,K] ~ normal(0, 1/psi);
    }  
  }

  for(i in 1:n){
    F[i] ~ multi_normal(mu_f,I);
  }

  for( i in 1:n) {
    for (j in 1:p) {
      Y[i,j] ~ normal(Lambda*F[i], sigma^2);
    }  
  }
}