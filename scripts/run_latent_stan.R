setwd("~/git.repos/BSFG/scripts/")
library("MASS")
library("rstan")

test <- generate_correlated_data(10, 6, 5,
          non_zero_entries = c(0.8, 0.7, 0.6, 0.5, 0.2))
test

n <- 10
p <- 6
K <- 5

I <- diag(rep(1, K))
I

Y <- test$Y

test.data <-list(n = n, p = p, K = K, I = I, Y = Y)
test.data

fit <- stan(file = 'simple_sparse_model.stan', data = test.data, 
            iter = 1000, chains = 4)