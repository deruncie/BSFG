library("rstan")
rstan_options(auto_write = T)
options(mc.cores = 1)


sim_dir = 'Sim_FE4_1'
load(sprintf('%s/Sim_data.RData',sim_dir))

# set priors
sim_data$K = 4
sim_data$nu = 3
sim_data$nu_B = 1
sim_data$alpha_B = 2.1
sim_data$beta_B = 1/5
sim_data$alpha1 = 2.1
sim_data$beta1 = 1/5
sim_data$alpha2 = 2.1
sim_data$beta2 = 1
sim_data$sigma_scale = 2.5
sim_data$F_vars_beta = 2

# now calculate Q matrix:
Z = sim_data$Z1
A = sim_data$A1
svd_ZAZ = svd(Z %*% A %*% t(Z))
sim_data$Q = svd_ZAZ$u
sim_data$d = svd_ZAZ$d
sim_data$d[sim_data$d< 1e-13] = 0

# if no X, add fake X
if(is.null(sim_data$X)){
	sim_data$X = matrix(1,nr=sim_data$n,nc = 1)
	sim_data$b = 1
}
# if no Z2, add fake Z2
if(is.null(sim_data$Z2)){
	sim_data$Z2 = matrix(0,nr=sim_data$n,nc = 0)
  sim_data$r2 = 0
  sim_data$A2_chol = diag(1,sim_data$r2)
}

Nitt = 100
warmup = 50
chains = 1

Full_model <- stan(file = 'Full_model.stan', chains = 0)

Full_model_fit <- sampling(object = get_stanmodel(Full_model), data = sim_data, 
            iter = Nitt, warmup = warmup,chains = chains, verbose = TRUE, refresh = 10
            ,control = list(
              # adapt_delta = 0.5,
              # max_treedepth = 6
              )
            ,pars = c("Lambda","QTF","F_vars","F_h2","E_h2","sigma2_a","sigma2_e","inv_tau","Y_hat","G","R","B","B_F","B_resid","mu","B_scale"), include = T
            )
fit = Full_model_fit


chi2_model <- stan(file = 'Full_model_repChi2.stan', chains = 0)

chi2_model_fit <- sampling(object = get_stanmodel(chi2_model), data = sim_data, 
            iter = Nitt, warmup = warmup,chains = chains, verbose = TRUE, refresh = 10
            ,control = list(
              # adapt_delta = 0.5,
              # max_treedepth = 6
              )
            ,pars = c("Lambda","QTF","F_vars","F_h2","E_h2","sigma2_a","sigma2_e","inv_tau","Y_hat","G","R","B","B_F","B_resid","mu","B_scale"), include = T
            )
fit = chi2_model_fit

simpleB_model <- stan(file = 'Full_model_simpleB.stan', chains = 0)

simpleB_model_fit <- sampling(object = get_stanmodel(simpleB_model), data = sim_data, 
            iter = Nitt, warmup = warmup,chains = chains, verbose = TRUE, refresh = 10
            ,control = list(
              # adapt_delta = 0.5,
              # max_treedepth = 6
              )
            ,pars = c("Lambda","QTF","F_vars","F_h2","E_h2","sigma2_a","sigma2_e","inv_tau","Y_hat","G","R","B","B_F","B_resid","mu","B_scale"), include = T
            )
fit = simpleB_model_fit

a = get_sampler_params(fit, inc_warmup = T)[[1]]
summary(do.call(rbind, args = list(a)), digits = 2)
summary(do.call(rbind, args = get_sampler_params(fit, inc_warmup = F)), digits = 2)

# compare fitted values to true values
y_hat = get_posterior_mean(fit,pars='Y_hat')
plot(t(sim_data$Y),y_hat);abline(0,1)

# compare fitted Lambda to true Lambda
Lambda_est = matrix(get_posterior_mean(fit,pars='Lambda'),sim_data$p,sim_data$K,byrow=F)
Lambda_act = sim_data$Lambda
corL = cor(Lambda_est,Lambda_act)

# re-order Lambda_act to correspond to Lambda_est
matched_cols = 1:4#apply(abs(corL),1,function(x) which(x==max(x)))
col_signs = sapply(1:length(matched_cols),function(i) sign(corL[i,matched_cols[i]]))
Lambda_act = sweep(Lambda_act[,matched_cols],2,col_signs,'*')
corL = corL[,matched_cols]

image(corL^2,zlim=c(0,1))
corL
apply(abs(corL),1,max)

plot(Lambda_act,Lambda_est,col = matrix(1:sim_data$K,nr=sim_data$p,nc = sim_data$K,byrow=T));abline(0,1)
plot(sim_data$Lambda %*% t(sim_data$Lambda),Lambda_est  %*% t(Lambda_est));abline(0,1)

G = matrix(get_posterior_mean(fit,pars='G'),nc=sim_data$p)
R = matrix(get_posterior_mean(fit,pars='R'),nc=sim_data$p)
B = matrix(get_posterior_mean(fit,pars='B'),nc=sim_data$p,byrow=T)
B_F = matrix(get_posterior_mean(fit,pars='B_F'),nc=sim_data$K_GT,byrow=T)
F_vars = matrix(get_posterior_mean(fit,pars='F_vars'),nc=3,byrow=T)

# plot(G,sim_data$G);abline(0,1)
# plot(R,sim_data$R);abline(0,1)
plot(B_F %*% t(Lambda_est),sim_data$B_F %*% t(sim_data$Lambda));abline(0,1)
plot(B,sim_data$B + sim_data$B_F %*% t(sim_data$Lambda));abline(0,1)
plot(B_F,sim_data$B_F);abline(0,1)
plot(get_posterior_mean(fit,pars='E_h2')[,1],sim_data$E_h2);abline(0,1)
plot(get_posterior_mean(fit,pars='mu')[,1],sim_data$mu);abline(0,1)

cor(c(B),c(sim_data$B + sim_data$B_F %*% t(sim_data$Lambda)))


e = monitor(fit,print=F)
e[grep('B_F',rownames(e)),]
e[grep('B[',rownames(e),fixed=T),]
e[grep('F_h2',rownames(e)),]
image(matrix(e[grep('Lambda[',rownames(e),fixed=T),'n_eff'],nc=sim_data$K,byrow=T))

print(fit)
plot(fit)
rstan::traceplot(fit, pars = "inv_tau", inc_warmup = T)
rstan::traceplot(fit, pars = "Lambda", inc_warmup = T)
rstan::traceplot(fit, pars = "F_h2", inc_warmup = T)
rstan::traceplot(fit, pars = "E_h2", inc_warmup = FALSE)
rstan::traceplot(fit, pars = "B_F", inc_warmup = FALSE)
rstan::traceplot(fit, pars = "B", inc_warmup = FALSE)
rstan::traceplot(fit, pars = "F_vars", inc_warmup = T)
plot(fit, pars = "inv_tau")
plot(fit, pars = "F_h2")
plot(fit, pars = "E_h2")
plot(fit, pars = c("B_resid",'B_F'))
plot(fit, pars = c('B'))
plot(fit, pars = c('B_F'))
plot(fit, pars = c('mu'))
plot(fit, pars = c('B_scale'))
plot(fit, pars = c('F_vars'))


e = extract(fit,inc_warmup=TRUE)
pairs(fit,pars='inv_tau')
pairs(fit,pars='F_h2',inc_warmup=T)

plot(fit)
rstan::traceplot(fit)

# decrease the step size
# delta1 min be 1
# sigma ~ cauchy(0, 1);

plot(fit_v3, pars = c('B'))
plot(BSFG_fit, pars = c('B'))

e_v3 = monitor(fit_v3,print=F)
e_v3[grep('B[',rownames(e_v3),fixed=T),'n_eff']
summary(e_v3[grep('B[',rownames(e_v3),fixed=T),'75%']-e_v3[grep('B[',rownames(e_v3),fixed=T),'25%'])
e_BSFG = monitor(BSFG_fit,print=F)
e_BSFG[grep('B[',rownames(e_BSFG),fixed=T),'n_eff']
summary(e_BSFG[grep('B[',rownames(e_BSFG),fixed=T),'75%']-e_BSFG[grep('B[',rownames(e_BSFG),fixed=T),'25%'])


library(ggplot2)
library(reshape)
e = extract(fit,pars='Lambda')[[1]]
means = apply(e,c(1,3),function(x) mean(x))
means = melt(means)
p = ggplot(means,aes(x=iterations,y=value,group=Var.2)) + geom_line(aes(color=Var.2));p