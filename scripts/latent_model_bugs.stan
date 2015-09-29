// https://groups.google.com/forum/#!topic/stan-users/qGAUsyJWyMA

data{
int N;          // sample size
int P;          // number of variables
int D;          // number of dimensions
vector[P] X[N]; // data matrix of order [N,P]
}

parameters{
vector[P] b;            // intercepts
vector<lower=0>[P] lam; // factor loadings
vector[D] FS[N];        // factor scores, matrix of order [N,D]
corr_matrix[D] Rho;     // correlation matrix between factors
vector<lower=0,upper=100>[P] var_p; // variance for each variable
}

transformed parameters{
vector[D] M;
vector<lower=0, upper=1000>[D] Sd_d; // sd of factors
vector[P] mu[N];
matrix[D,D] Ld;

for (m in 1:D) {
M[m] <- 0;
Sd_d[m] <- 1;}

Ld <- diag_matrix(Sd_d) * cholesky_decompose(Rho);

for(i in 1:N){
    mu[i,1] <- b[1] + lam[1]*FS[i,1];
    mu[i,2] <- b[2] + lam[2]*FS[i,1];
    mu[i,3] <- b[3] + lam[3]*FS[i,1];
    mu[i,4] <- b[4] + lam[4]*FS[i,2];
    mu[i,5] <- b[5] + lam[5]*FS[i,2];
    mu[i,6] <- b[6] + lam[6]*FS[i,2];
    mu[i,7] <- b[7] + lam[7]*FS[i,3];
    mu[i,8] <- b[8] + lam[8]*FS[i,3];
    mu[i,9] <- b[9] + lam[9]*FS[i,3];
    }
}

model{

b ~ normal(0, 100);
lam ~ normal(0.4,var_p);

for(i in 1:N){   
    X[i] ~ normal(mu[i],var_p);
    FS[i] ~ multi_normal_cholesky(M, Ld);
    }

}