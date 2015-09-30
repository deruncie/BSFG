setwd("~/git.repos/BSFG/scripts/")
library("MASS")
library("rstan")
library("parallel")
rstan_options(auto_write = T)
options(mc.cores = 4)


n <- 10
p <- 60
K <- 5
test <- generate_correlated_data(n,p,K,
          non_zero_entries = c(0.3, 0.2, 0.15, 0.10, 0.1),
          ideosyncratic_variance_parms = c(20,1))
test


I <- diag(rep(1, K))
I


Y <- test$Y
nu <- 3.0
alpha1 <- 2
alpha2 <- 2


test.data <-list(n = n, p = p, K = K, I = I, Y = Y, nu = nu,
	              alpha1 = alpha1, alpha2 = alpha2)
test.data


test_model <- stan(file = 'simple_sparse_model.stan', chains = 0)

fit <- sampling(object = get_stanmodel(test_model), data = test.data, 
            iter = 1000, chains = 1, verbose = TRUE, refresh = 10,
            control = list(adapt_delta = 0.95),
            pars = c("Lambda","Y_hat", "tau"), include = FALSE)

y_hat = get_posterior_mean(fit,pars='Y_hat')

plot(t(Y),y_hat)
l_hat = array(get_posterior_mean(fit,pars='Lambda'),dim = dim(t(test$Lambda)))
l_hat = t(l_hat)
cor(l_hat,test$Lambda)
 summary(do.call(rbind, args = get_sampler_params(fit, inc_warmup = FALSE)), digits = 2)


plot(test$Lambda,l_hat);abline(0,1)

plot(test$Lambda %*% t(test$Lambda),l_hat  %*% t(l_hat));abline(0,1)

print(fit)
plot(fit)
traceplot(fit, pars = "tau", inc_warmup = FALSE)
plot(fit, pars = "tau")
# decrease the step size
# delta1 min be 1
# sigma ~ cauchy(0, 1);