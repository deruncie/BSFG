setwd("~/git.repos/BSFG/scripts/")
library("MASS")
library("rstan")
library("parallel")
rstan_options(auto_write = T)
options(mc.cores = 4)


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


test_model <- stan(file = 'simple_sparse_model.stan', chains = 0)

fit <- sampling(object = get_stanmodel(test_model), data = test.data, 
            iter = 1000, chains = 1, verbose = TRUE, refresh = 10, 
            pars = c("Lambda","Y_hat"), include = FALSE)

print(fit)
plot(fit)
traceplot(fit)