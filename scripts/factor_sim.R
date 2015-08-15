# Modified slightly from example by Rick Farouni
# http://rfarouni.github.io/assets/projects/BayesianFactorAnalysis/BayesianFactorAnalysis.html
# 
library("MASS")
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



