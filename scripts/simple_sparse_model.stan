// simple sparse factor model
// latent portion
// Cody Markelz -- markelz@gmail.com

functions{
  matrix sqrt_mat(matrix X){
    matrix[rows(X),cols(X)] res;
    for (i in 1:rows(X)){
      for(j in 1:cols(X)){
        res[i,j] <- sqrt(X[i,j]);
      }
    }
    return(res);
  }

  vector cum_prod(vector X){
    vector[rows(X)] res;
    res[1] <- X[1];
    for (i in 2:rows(X)){
      res[i] <- res[i-1]*X[i];
    }
    return(res);
  }
}

data {
  int<lower=1> n;                // number of individuals 1 ...n
  int<lower=1> p;                // number of traits 1 ... p
  int<lower=1> K;                // number of latent factors 
  matrix[n,p]  Y;                // data matrix of order [n,p]
  matrix[K,K]  I;
  real          nu;
  real         alpha1;
  real         alpha2;

}

transformed data {
  vector[K] mu_f;
  mu_f <- rep_vector(0.0, K);   // mu_f zero vector
}

parameters{
  vector[K] F[n];
  vector<lower=0>[p]   sigma;
  matrix<lower=0>[p,K] psi;
  matrix[p,K]          Lambda_std;
  vector<lower=0>[1]    delta;
  vector<lower=0>[K-1] delta1;

  //real            lam;
}

transformed parameters{
  vector[K] tau;
  matrix[p,K] Lambda;
  tau <- cum_prod(append_row(delta,delta1));
  Lambda <- Lambda_std ./ sqrt_mat(diag_post_multiply(psi,tau));

}

model {
  vector[p] mu;
  for(i in 1:p){
    row(Lambda_std,i) ~ normal(0,1);
    row(psi,i) ~ gamma(nu/2.0, nu/2.0);
  }

  delta ~ gamma(alpha1, 1);
  delta1 ~ gamma(alpha2,1);

  sigma ~ cauchy(0, 1);

  for(i in 1:n){
    F[i] ~ multi_normal(mu_f,I);
  }

  mu <- Lambda*F[1];
  //print(mu);
  for( i in 1:n) {
    //mu <- Lambda*F[i];
    row(Y,i) ~ normal(Lambda*F[i], sigma);
  }
}

generated quantities {
  vector[p] Y_hat[n];

  for( i in 1:n) {
    //mu <- Lambda*F[i];
    Y_hat[i] <- Lambda*F[i];
  }
}
