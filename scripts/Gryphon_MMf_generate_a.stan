functions {

  vector log_lik_normal(vector y, vector mu, vector sd){
    vector[rows(y)] res;
    for(i in 1:rows(y)){
      res[i] <- normal_log(y[i],mu[i],sd[i]);
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
	vector sample_a_rng(vector y, matrix UtZt, vector d, matrix LiU){
		int r;
		vector[rows(d)] mu;
		vector[rows(d)] a_tilde;
		vector[rows(d)] a;
		r <- rows(d);
		mu <- diag_pre_multiply(d,UtZt) * y;
		for(i in 1:r){
			a_tilde[i] <- normal_rng(mu[i],sqrt(d[i]));
		}
		a <- LiU*a_tilde;
		return(a);
	}
}
data {
	int<lower=1> N;
	int r1;
	vector[N] Y;
	matrix[N,r1] Z1;
	int<lower=0> p_X;
	matrix[N,p_X] X;
	matrix[N,N] U;
	vector[N] D;
	matrix[r1,N] UtZt;
	vector[r1] U_D;
	matrix[r1,r1] LiU;
}
transformed data{
	vector[N] UtY;
	matrix[N,p_X] UtX;

	 UtY <- U'*Y; //'
	 UtX <- U'*X; //'
}
parameters {
	vector[p_X] beta;
	real<lower=0,upper=0.99> h2;
	real<lower=0, upper=pi()/2> sigma_p_unif;
}
transformed parameters {
	real sigma_p;

	sigma_p <- 5 * tan(sigma_p_unif);
}
model {
	real sigma;
	real sigma_Z1;
	
	beta ~ normal(0,1000);

	sigma <- sqrt((1-h2)*sigma_p^2);
	sigma_Z1 <- sqrt(h2*sigma_p^2);

	UtY ~ normal(UtX*beta,sqrt_vec(sigma_Z1^2*D + sigma^2));
}
generated quantities{
	real sigma;
	real sigma_Z1;
	real sigma2;
	real sigma_Z12;
	real sigma_p2;
	vector[r1] a;

	sigma <- sqrt((1-h2)*sigma_p^2);
	sigma_Z1 <- sqrt(h2*sigma_p^2);
	
	sigma2 <- sigma^2;
	sigma_Z12 <- sigma_Z1^2;
	sigma_p2 <- sigma_p^2;

	a <- sample_a_rng(Y-X*beta,UtZt,1.0/sigma_Z12 * U_D + 1.0/sigma2,LiU);
}