// genetic sparse factor model
// Based on Runcie and Mukherjee 2013
// Cody Markelz -- markelz@gmail.com
// Daniel Runcie --deruncie@ucdavis.edu
// v2c - adding fixed effects for F
// v3 - altering fixed effect structure for F vs Y
// 

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
  matrix[r2,r2] A2_chol;         // cholesky factor for A2 such that A2_chol * A2_chol' = A2
  real         nu_B;      // shrinkage for Bs
  real         alpha_B;
  real         beta_B;
  real         nu;        // shrinkage for Lambdas
  real         alpha1;
  real         beta1;
  real         alpha2;
  real         beta2;
  real         sigma_scale;
  matrix[n,n]  Q;
  vector[n]    d;
  int          F_vars_beta;
}

transformed data {
  matrix[n,n]   QT;
  vector[n]     QT_mu;
  vector[n]     QTY[p];
  matrix[n,b]   QTX;
  matrix[n,r2]  QTZ2;
  int           num_vars; // 2 + any_fixed_effects + any_extra_random_effects
  print("N: ",n," p: ",p," b: ",b," K: ",K);
  
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

  num_vars <- 2;
  if(b > 0) num_vars <- num_vars + 1;
  if(r2 > 0) num_vars <- num_vars + 1;
  print(num_vars);
}

parameters{
  row_vector[p] mu;       // global mean for each gene
  matrix[r2,p]  U2_std;   // location effects of random effect 2

  vector[b]             B_scale;        // scale factor for B coefficients
  matrix[b,p]           B_std;          // std_normal underlying fixed effect coefficients for Y
  matrix<lower=0>[b,p]  B_psi;          // precision of fixed effect coefficients for Y

  vector[n]             F_std[K];       // std_normal underlying F
  vector[b]             B_F_std[K];     // std_normal underlying fixed effect coefficients for F
  vector[r2]            U2_F_std[K];     // std_normal underlying random effect 2 coefficients for F
  vector[F_vars_beta]   F_vars_z[K,num_vars];  // components of re-parameterized simplex F_vars - a,r,b,a2

  vector<lower=0,upper=pi()/2>[p] sigma_unif[2+(r2>0)];  // variances for residuals: 1: A, 2: I, 3: A2 - re-parameterized cauchy

  row_vector[p]          Lambda_std[K]; // std_normal underlying Lambda mixture
  row_vector<lower=0>[p] psi[K];        // part of the precision mixture contributing to t-distribution on Lambda components
  vector<lower=0>[K]     delta;    // shrinkage components of tau
}

transformed parameters{
  vector[p]     sigma_U2;   // residual genetic sd
  vector[p]     sigma2_a;      // total variation (G+E) of each gene's residuals on the factors
  vector[p]     sigma2_e;      // total variation (G+E) of each gene's residuals on the factors
  row_vector[p] Lambda[K]; 
  vector[K]     tau;  
  vector[n]     QTF[K];     // rotated factor scores  
  matrix[b,p]   B_resid;    // fixed effect coefficients for Y
  vector[b]     B_F[K];     // fixed effect coefficients for F
  vector[r2]    U2_F[K];    // random effect 2 coefficients for F
  simplex[num_vars]    F_vars[K];   // variances of variance components of F (fixed, random, resid)

// sds - re-parmeterized cauchy
  sigma2_a <- pow_vec(sigma_scale * tan_vec(sigma_unif[1]),2.0);
  sigma2_e <- pow_vec(sigma_scale * tan_vec(sigma_unif[2]),2.0);
  if(r2 > 0) {
    sigma_U2 <- sigma_scale * tan_vec(sigma_unif[3]);
  }

// Fixed effect residuals
  if(b > 0){
    B_resid <- B_std ./ sqrt_mat(diag_pre_multiply(B_scale,B_psi));
  }

// Factors
  if(K > 0){
// Lambda - factor loadings
    tau <- cum_prod(delta);
    for(k in 1:K){
      Lambda[k] <- Lambda_std[k] ./ sqrt_rvec(psi[k] * tau[k]);
    }

  // Factor scores
    // calculate simplex F_vars.
    for(k in 1:K){
      for(l in 1:num_vars){
        F_vars[k][l] <- dot_self(F_vars_z[k,l]);
      }
      F_vars[k] <- F_vars[k] / sum(F_vars[k]);
    }

    // Build up F column by column
    for(k in 1:K){
      QTF[k] <- F_std[k] .* sqrt_vec(F_vars[k][1]*d + F_vars[k][2]); // standard deviation of QTF from random + resid
    }
    if(b > 0){
      for(k in 1:K){
        B_F[k] <- B_F_std[k] * sqrt(F_vars[k][3]);
        QTF[k] <- QTF[k] + QTX * B_F[k];
      }
    }
    if(r2 > 0){
      int var_index;
      var_index <- 3 + (b>0);
      for(k in 1:K){
        U2_F[k] <- A2_chol * U2_F_std[k] * sqrt(F_vars[k][var_index]);
        QTF[k] <- QTF[k] + QTZ2 * U2_F[k];
      }
    }
  }
}

model {
  vector[n] sd_E[p];  // residual standard deviations - including main random effect

// note: mu has flat priors
  mu ~ normal(0,1000);

// Random effect 2
  if(r2 > 0) {
    for(i in 1:n){    
      U2_std[i] ~ normal(0,1);
    }
    if(K > 0){
      for(k in 1:K){
        U2_F_std[k] ~ normal(0,1);        
      }
    }
  }

// components of B_resid
  // shrinkage
  if(b > 0){
    B_scale ~ gamma(alpha_B,beta_B);
    // coefficients
    for(x in 1:b){
      B_std[x] ~ normal(0,1);
      B_psi[x] ~ gamma(nu_B/2.0,nu_B/2.0);
    }
  }

//Factors
  if(K > 0){
  // components of Lambda
    // column-wise shrinkage
    delta[1] ~ gamma(alpha1, beta1); // global shrinkage
    if(K>1){
      segment(delta,2,K-1) ~ gamma(alpha2, beta2); // additional shrinkage
    }
    // coefficients
    for(k in 1:K){
      Lambda_std[k] ~ normal(0,1);
      psi[k] ~ gamma(nu/2.0, nu/2.0);
    }
    
  // components of F
    for(k in 1:K){
      F_std[k] ~ normal(0,1);
      B_F_std[k] ~ normal(0,1);
    }   
    // individual components of F_vars simplex
    for(k in 1:K){
      for(l in 1:num_vars) {
        F_vars_z[k,l] ~ normal(0,1);
      }
    }
  }
  
  //calculate rotated residual standard deviations - including main random effect
  for(j in 1:p){
    sd_E[j] <- sqrt_vec(sigma2_a[j] * d + sigma2_e[j]);
  }

// likelihood
  {
    matrix[n,p] QTY_mean;
    QTY_mean <- QT_mu * mu;
    if(b > 0){
      QTY_mean <- QTY_mean + QTX * B_resid;
    }
    // add in random effect 2
    if(r2 > 0){
      QTY_mean <- QTY_mean + QTZ2 * A2_chol * diag_post_multiply(U2_std,sigma_U2);
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
  vector[K] F_h2;
  matrix[p,p] G;
  matrix[p,p] R;
  matrix[b,p] B;
  
  inv_tau <- rep_vector(1.0,K) ./ tau;
  
  Y_hat <- rep_vector(1.0,n) * mu;
  G <- diag_matrix(sigma2_a);
  R <- diag_matrix(sigma2_e);
  B <- B_resid;

  if(b > 0){
    Y_hat <- Y_hat + X * B_resid;
  }

  if(r2 > 0){
    Y_hat <- Y_hat + Z2 * A2_chol * diag_post_multiply(U2_std,sigma_U2);
  }

  if(K > 0) {
    for(k in 1:K){
      Y_hat <- Y_hat + Q * QTF[k] * Lambda[k];
    }
    
    for(k in 1:K){
        F_h2[k] <- F_vars[k][1] / (F_vars[k][1]+F_vars[k][2]);
    }

    for(k in 1:K){
      G <- G + F_h2[k] * Lambda[k]' * Lambda[k]; //'
      R <- R + (1.0 - F_h2[k]) * Lambda[k]' * Lambda[k]; //'
    }

    for(k in 1:K){
      B <- B + B_F[k] * Lambda[k];
    }
  }
  
}
