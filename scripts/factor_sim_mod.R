# Modified slightly from example by Rick Farouni
# http://rfarouni.github.io/assets/projects/BayesianFactorAnalysis/BayesianFactorAnalysis.html
# 
setwd("~/git.repos/BSFG/scripts/")
library("MASS")
library("rstan")
library("parallel")
set.seed(42)
D <- 3 # number of factors
P <- 10 # number of predictors
N <- 300 # number of individuals

mu_theta <- rep(0, D)
mu_theta

mu_epsilon <- rep(0, P)

Phi <- diag(rep(1, D))
Phi
#      [,1] [,2] [,3]
# [1,]    1    0    0
# [2,]    0    1    0
# [3,]    0    0    1

Psi <- diag(c(0.2079, 0.19, 0.1525, 0.20, 0.36, 0.1875, 0.1875, 1.00, 0.27, 0.27))
Psi 

# individual loadings for factors
l1 <- c(0.99, 0.00, 0.25, 0.00, 0.80, 0.00, 0.50, 0.00, 0.00, 0.00)
l2 <- c(0.00, 0.90, 0.25, 0.40, 0.00, 0.50, 0.00, 0.00, -0.30, -0.30)
l3<-  c(0.00, 0.00, 0.85, 0.80, 0.00, 0.75, 0.75, 0.00, 0.80, 0.80)

L <- cbind(l1, l2, l3) # loading matrix
L

?mvrnorm
theta <- mvrnorm(N, mu_theta, Phi)
theta

epsilon <- mvrnorm(N, mu_epsilon, Psi)
dim(epsilon)

Y <- theta%*%t(L) + epsilon # observations
dim(Y)

M  <- D*(P-D)+ D*(D-1)/2
M

fa.data <-list(P=P,N=N,Y=Y,D=D)
fa.data2 <-list(P=P,N=N,Y=Y)

# a function to generate intial values that are slightly jittered for each chain.
init_fun = function() {
  init.values<-list(L_t=rep(0,24)+runif(1,-.1,.1),
                    L_d=rep(.5,D)+runif(1,-.1,.1),
                    psi=rep(.2,P)+runif(1,-.1,.1),
                    sigma_psi=0.15+runif(1,-.1,.1),
                    mu_psi=0.2++runif(1,-.1,.1),
                    sigma_lt=0.5+runif(1,-.1,.1),
                    mu_lt=0.0+runif(1,-.1,.1))
  return(init.values); 
} 

#compile the model
fa.model<- stan("latent_model_mod2.stan", 
                  data = fa.data,
                  chains =1, 
                  pars=c("L","psi","sigma_psi","mu_psi","sigma_lt","mu_lt"))

# run 4 chain in parallel and save samples onto file
Nchains <- 1
Niter <- 300
t_start <- proc.time()[3]
sflist <- mclapply(1:Nchains, mc.cores = Nchains, 
                     function(i) stan(fit = fa.model, 
                                      data =fa.data, 
                                      pars= c("L","psi","sigma_psi","mu_psi","sigma_lt","mu_lt"), 
                                      seed = 42,
                                      iter=Niter,
                                      init=init_fun,
                                     #diagnostic_file = paste("diagfile",i,".csv",sep = ""),
                                      sample_file = paste("sampfile",i,".csv",sep = ""),
                                      chains = 1, chain_id = i, 
                                      refresh = -1))

t_end <- proc.time()[3]
t_elapsed <- t_end - t_start
 
fa.fit<- sflist2stanfit(sflist) 
print(fa.fit,probs = c(0.5))






