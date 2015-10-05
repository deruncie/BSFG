// model from Runcie and Mukherjee 2013
// Cody Markelz -- markelz@gmail.com
// Daniel Runcie --deruncie@ucdavis.edu

functions{
  vector tan_vec(vector X){
    vector[rows(X)] res;
    for(i in 1:rows(X)){
        res[i] <- tan(X[i]);
    }
    return(res);
  }
  vector pow_vec(vector X,real n){
    vector[rows(X)] res;
    for(i in 1:rows(X)){
        res[i] <- pow(X[i],n);
    }
    return(res);
  }
  vector sqrt_vec(vector X){
    vector[rows(X)] res;
    for(i in 1:rows(X)){
        res[i] <- sqrt(X[i]);
    }
    return(res);
  }
  row_vector sqrt_rvec(row_vector X){
    row_vector[cols(X)] res;
    for(i in 1:cols(X)){
        res[i] <- sqrt(X[i]);
    }
    return(res);
  }
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
  int<lower=0> K;                // number of latent factors 
  int<lower=1> b;                // number of fixed effects 
  int<lower=0> r2;               // number of levels of random effect 2 (not in Q)
  matrix[n,p]  Y;                // data matrix of order [n,p]
  matrix[n,b]  X;                // fixed effect design matrix
  matrix[n,r2] Z2;               // random effect 2 design matrix
  real         nu;        // shrinkage for Lambdas
  real         alpha1;
  real         beta1;
  real         alpha2;
  real         beta2;
  real         sigma_scale;
  matrix[n,n]  Q;
  vector[n]    d;
}

transformed data {
  matrix[n,n]   QT;
  vector[n]     QTY[p];
  vector[n]     QT_mu;
  matrix[n,b]   QTX;
  matrix[n,r2]  QTZ2;
  print("N: ",n," p: ",p," b: ",b," r2: ",r2," K: ",K);
  
  QT <- Q'; //'

  // rotate Y by QT. Extract rows for quicker axis
  {
    for(j in 1:p){
      QTY[j] <- QT * col(Y,j);
    }
  }
  
  // rotate mean mu;
  QT_mu <- QT * rep_vector(1.0,n);

  // rotate X by QT
  QTX <- QT*X;

  // rotate Z2 by QT
  if(r2 > 0){
    QTZ2 <- QT*Z2;
  }
}

parameters{
  row_vector[p] mu;     // global mean for each gene
  matrix[b,p]   B;      // fixed effect coefficients for Y
  matrix[r2,p]  U2_std;     // location effects of random effect 2

  vector[n]                    F_std[K];       // std_normal underlying F
  vector<lower=0,upper=.95>[K] F_h2;            // residual heritability

  vector<lower=0,upper=pi()/2>[p] sigma_U2_unif;      // sd 2nd random effect - re-parameterized cauchy
  vector<lower=0,upper=pi()/2>[p] sigma_a_unif;      // residual genetic sd - re-parameterized cauchy
  vector<lower=0,upper=pi()/2>[p] sigma_e_unif;      // residual residual sd - re-parameterized cauchy

  row_vector[p]          Lambda_std[K]; // std_normal underlying Lambda mixture
  row_vector<lower=0>[p] psi[K];        // part of the precision mixture contributing to t-distribution on Lambda components
  vector<lower=0>[1]     delta_1;    // first shrinkage component of tau
  vector<lower=0>[K-1]   delta_K;    // remaining shrinkage components of tau

}

transformed parameters{
  vector[p]     sigma_U2;      // residual genetic sd
  vector[p]     sigma2_a;      // residual genetic variance
  vector[p]     sigma2_e;      // residual residual variance
  row_vector[p] Lambda[K]; 
  vector[K]     tau;  
  vector[n]     QTF[K];     // rotated factor scores  

// total residual variance - re-parmeterized cauchy
  sigma_U2 <- sigma_scale * tan_vec(sigma_U2_unif);
  sigma2_a <- pow_vec(sigma_scale * tan_vec(sigma_a_unif),2.0);
  sigma2_e <- pow_vec(sigma_scale * tan_vec(sigma_e_unif),2.0);

// Factors
  if(K > 0){
// Lambda - factor loadings
    tau <- cum_prod(append_row(delta_1,delta_K));
    for(k in 1:K){
      Lambda[k] <- Lambda_std[k] ./ sqrt_rvec(psi[k] * tau[k]);
    }

  // Factor scores
    // Build up F column by column
    for(k in 1:K){
      vector[n] sd_F;  // standard deviation of QTF from random + resid
      sd_F <- sqrt_vec(F_h2[k]*d + (1.0-F_h2[k]));
      QTF[k] <- sd_F .* F_std[k];
    }
  }
}

model {
  vector[n] sd_E[p];  // residual standard deviations - including main random effect

// note: mu, B have flat priors

// Random effect 2
  if(r2 > 0) {
    for(i in 1:n){    
      U2_std[i] ~ normal(0,1);
    }
  }

//Factors
  if(K > 0){
  // components of Lambda
    // column-wise shrinkage
    delta_1 ~ gamma(alpha1, beta1); // global shrinkage
    delta_K ~ gamma(alpha2, beta2); // additional shrinkage
    // coefficients
    for(k in 1:K){
      Lambda_std[k] ~ normal(0,1);
      psi[k] ~ gamma(nu/2.0, nu/2.0);
    }
    
  // components of F
    for(k in 1:K){
      F_std[k] ~ normal(0,1);
    }   
    // F_h2 has flat priors
  }
  
  //calculate rotated residual standard deviations - including main random effect
  // E_h2 has flat priors
  for(j in 1:p){
    sd_E[j] <- sqrt_vec(sigma2_a[j] * d + sigma2_e[j]);
  }

// likelihood
  {
    matrix[n,p] QTY_mean;
    QTY_mean <- QT_mu * mu;
    if(b > 0){
      QTY_mean <- QTY_mean + QTX * B;
    }
    // add in random effect 2
    if(r2 > 0){
      QTY_mean <- QTY_mean + QTZ2 * diag_post_multiply(U2_std,sigma_U2);
    }
    // add in factors
    if(K > 0){
      for(k in 1:K){
        QTY_mean <- QTY_mean + QTF[k] * Lambda[k]; 
      }
    }
    for(j in 1:p) {
      QTY[j] ~ normal(col(QTY_mean,j),sd_E[j]);
    }
  }
}

generated quantities {
  matrix[n,p] Y_hat;
  vector[K] inv_tau;
  vector[p] E_h2;
  matrix[p,p] G;
  matrix[p,p] R;
  
  inv_tau <- rep_vector(1.0,K) ./ tau;
  E_h2 <- sigma2_a ./ (sigma2_a + sigma2_e);
  
  Y_hat <- rep_vector(1.0,n) * mu;
  G <- diag_matrix(sigma2_a);
  R <- diag_matrix(sigma2_e);

  if(b > 0){
    Y_hat <- Y_hat + X * B_resid;
  }

  if(r2 > 0){
    Y_hat <- Y_hat + Z2 * diag_post_multiply(U2_std,sigma_U2);
  }

  if(K > 0) {
    for(k in 1:K){
      Y_hat <- Y_hat + Q * QTF[k] * Lambda[k];
    }
    
    for(k in 1:K){
      G <- G + F_h2[k] * Lambda[k]' * Lambda[k]; //'
      R <- R + (1.0 - F_h2[k]) * Lambda[k]' * Lambda[k]; //'
    }
  }
  
}
