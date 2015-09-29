// simple sparse factor model
// latent portion
// Cody Markelz -- markelz@gmail.com

data {


}

transformed data {


}

parameters{


}

transformed parameters {

}

model {
  psi ~ gamma(3/2, 3/2);
  sigma ~ cauchy(0, 1);

  f[i] ~ normal(0,I)
  y[i] ~ normal(mu, sig^2)

}