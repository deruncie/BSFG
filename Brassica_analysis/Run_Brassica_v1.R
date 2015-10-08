library("rstan",lib.loc='~/R/x86_64-pc-linux-gnu-library/3.2/')
rstan_options(auto_write = T)
options(mc.cores = 1)

in_data = read.csv('rlog2_shade_brassica_shade.csv')

data = data.frame(ID = colnames(in_data)[-c(1)])
data$Genotype = sapply(as.character(data$ID),function(x) paste(strsplit(x,'_')[[1]][1:2],collapse='_'))
data$TRT = sapply(as.character(data$ID),function(x) strsplit(x,'_')[[1]][3])

Y = t(in_data[,-1])
Y = sweep(Y,1,rowMeans(Y),'-')
Y = sweep(Y,1,apply(Y,1,sd),'/')
# data = data[1:400,]
# Y = Y[1:400,1:50]

Bra_data = list()
Bra_data$Y = Y
Bra_data$Z1 = model.matrix(~Genotype+0,data)
Bra_data$X = matrix(model.matrix(~TRT,data)[,-1],nc=1)
Bra_data$A1 = diag(1,ncol(Bra_data$Z1))
Bra_data$n = nrow(Bra_data$Y)
Bra_data$p = ncol(Bra_data$Y)
Bra_data$b = ncol(Bra_data$X)

# set priors
Bra_data$K = 10
Bra_data$nu = 5
Bra_data$nu_B = 5
Bra_data$alpha_B = 2.1
Bra_data$beta_B = 1/5
Bra_data$alpha1 = 2.1
Bra_data$beta1 = 1/5
Bra_data$alpha2 = 2.1
Bra_data$beta2 = 1
Bra_data$sigma_scale = 2.5
Bra_data$F_vars_beta = 2

# now calculate Q matrix:
Z = Bra_data$Z1
A = Bra_data$A1
svd_ZAZ = svd(Z %*% A %*% t(Z))
Bra_data$Q = svd_ZAZ$u
Bra_data$d = svd_ZAZ$d
Bra_data$d[Bra_data$d< 1e-13] = 0

# if no X, add fake X
if(is.null(Bra_data$X)){
	Bra_data$X = matrix(1,nr=Bra_data$n,nc = 1)
	Bra_data$b = 1
}
# if no Z2, add fake Z2
if(is.null(Bra_data$Z2)){
	Bra_data$Z2 = matrix(0,nr=Bra_data$n,nc = 0)
	Bra_data$r2 = 0
  Bra_data$A2_chol = matrix(0,0,0)
}

Nitt = 300
warmup = 200
chains = 1

Full_model <- stan(file = '../scripts/Full_model.stan', chains = 0)

Full_model_fit <- sampling(object = get_stanmodel(Full_model), data = Bra_data, 
            iter = Nitt, warmup = warmup,chains = chains, verbose = TRUE, refresh = 10
            ,control = list(
              # adapt_delta = 0.5,
              # max_treedepth = 12
              )
            ,pars = c("Lambda","QTF","F_vars","F_h2","E_h2","sigma2_a","sigma2_e","inv_tau","Y_hat","G","R","B","B_F","mu","B_scale"), include = T
            )
fit = Full_model_fit
save(fit,file = sprintf('Bra_data_fit_K_%d.RData',sim_data$K))

