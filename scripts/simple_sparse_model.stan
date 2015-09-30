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
}

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
  vector<lower=0>[p]   sigma;
  matrix<lower=0>[p,K] psi;
  matrix[p,K]          Lambda_std;
  //real            lam;
}

transformed parameters{
//  diag_matrix(rep_vector(1, K)) I; // factor identity matrix
  matrix[p,K] Lambda;
  //Lambda <- diag_pre_multiply(Lambda_std * sqrt(1/ps);
  Lambda <- Lambda_std ./ sqrt_mat(psi);
  //print(row(Lambda,1));
//Lambda <- operator./(Lambda_std,sqrt_mat(psi))

}

model {
  vector[p] mu;
  for(i in 1:p){
    row(Lambda_std,i) ~ normal(0,1);
    row(psi,i) ~ gamma(3.0/2.0, 3.0/2.0);
  }
  # Lambda_std ~ normal(0,1);
  sigma ~ cauchy(0, 1);
  



  // Lambda ~ normal(0,sqrt(1/pdi));
  // for( j in 1:p) {
  //   for (k in 1:K) {
  //     Lambda[j,K] ~ normal(0, 1/psi[j,k]);
  //   }  
  // }

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
