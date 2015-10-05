library(MCMCglmm)
data(BTped)
Nped <- BTped[which(apply(BTped, 1, function(x) {any(x == "R187920" | x == "R187921")})), ]


Aped <- 2 * kinship2::kinship(Nped[, 1], Nped[,2], Nped[, 3])
Aphylo <- vcv.phylo(Nphylo, cor = T)

Data <- as.data.frame(read.table(file = "./gryphon.txt", header = TRUE))
names(Data)[1] <- "animal"
Data$animal <- as.factor(Data$animal)
Data$MOTHER <- as.factor(Data$MOTHER)
Data$BYEAR <- as.factor(Data$BYEAR)
Data$SEX <- as.factor(Data$SEX)
Data$BWT <- as.numeric(Data$BWT)
Data$TARSUS <- as.numeric(Data$TARSUS)
head(Data)
Ped <- as.data.frame(read.table(file = "./gryphonped.txt", header = TRUE))
for (x in 1:3) Ped[, x] <- as.factor(Ped[, x])
head(Ped)





library(rstan)
rstan_options(auto_write = T)
options(mc.cores = 4)


stan_model = stan('Gryphon_MMf.stan',chains=0)
stan_model = stan('Gryphon_MMf_newNormal.stan',chains=0)
stan_model = stan('Gryphon_MMf_generate_a.stan',chains=0)

A = 2 * kinship2::kinship(Ped[,1],Ped[,2],Ped[,3])

# a = matrix(rnorm(nrow(A)*10),nc=10);a=a %*% t(a) + diag(10,ncol(A))
# colnames(a) = rownames(a) = colnames(A)
# A = diag(1/sqrt(diag(a))) %*% a %*% diag(1/sqrt(diag(a)))
# colnames(A) = rownames(A) = colnames(a)

i = which(!is.na(Data$BWT))[1:300]
A_index = match(Data$animal[i],colnames(A))
A_sub = A[A_index,A_index]
X = model.matrix(~1,Data[i,])
Z1 = diag(1,ncol(A_sub))
sZAZ = svd(Z1 %*% A_sub %*% t(Z1))
sZ = svd(Z1 %*% t(Z1))
A_chol = t(chol(A_sub))
stan_data = list(
	N = length(i),
	Y = Data$BWT[i],
	p_X = ncol(X),
	X = X,
	U = sZAZ$u,
	D = sZAZ$d,
  r1 = ncol(Z1),
  UtZt = t(sZ$u) %*% t(Z1),
  U_D = sZ$d,
  LiU = solve(A_chol) %*% sZ$u
	)


Nchains <- 1
Niter <- 2000
Nwarm <- 50

stan_fit = sampling(  object = get_stanmodel(stan_model), 
    data =stan_data, 
    iter=Niter,
    warmup = Nwarm,
    init='random',
    chains = Nchains, verbose=TRUE,
    refresh = 10,
    pars=c('beta','sigma2','sigma_Z12','h2','sigma_p'),include=F #
    )

a = get_posterior_mean(stan_fit,pars='a')
get_adaptation_info(stan_fit)
plot(stan_fit)

library(lme4)
lme1 = lmer(BWT~1+())

print(stan_fit)
pairs(stan_fit,pars=c('h2','sigma_p'))
e = do.call(cbind,extract(stan_fit,pars=c('sigma','sigma_Z1','h2')))
h2 = e[,2]^2/(e[,2]^2+e[,1]^2)
plot(h2,e[,3])
abline(0,1)

# now using sparse Z1
stan_model2 = stan('Gryphon_MMg.stan',chains=0)
sparse_Z1 = as(Z1,'dgRMatrix')
stan_data2 = list(
  N = length(i),
  Y = Data$BWT[i],
  p_X = ncol(X),
  X = X,
  p_Z1 = ncol(Z1),
  non_zero_Z1 = sum(Z1 != 0),
  w_Z1 = sparse_Z1@x,
  v_Z1 = sparse_Z1@j+1,
  u_Z1 = sparse_Z1@p+1
)

stan_fit2 = sampling(  object = get_stanmodel(stan_model2), 
                      data =stan_data2, 
                      iter=Niter,
                      warmup = Nwarm,
                      init='random',
                      chains = Nchains, verbose=TRUE,
                      refresh = 10,
                      pars=c('beta','sigma2','sigma_Z12','h2','sigma_p') #
)

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
model1.1 <- MCMCglmm(BWT ~ 1, random = ~animal, pedigree = Ped, data = Data[i,], nitt = 6500, thin = 50, burnin = 1500, prior = prior1.1, verbose = FALSE,pr=T,pl=T)
summary(model1.1)

m = colMeans(model1.1$Sol)[paste0("animal.",colnames(A_sub))]
plot(m,a);abline(0,1)


sampler_params <- get_sampler_params(stan_fit, inc_warmup = T)
names(sampler_params) <- colnames(stan_fit)

lapply(sampler_params, FUN = function(x) {step <- table(x[,"stepsize__"]); names(step)[step == max(step)]})
lapply(sampler_params, FUN = function(x) table(x[,"treedepth__"]))
lapply(sampler_params, FUN = function(x) mean(x[,"accept_stat__"]))
lapply(sampler_params, FUN = function(x) table(x[,"n_divergent__"]))
lapply(sampler_params, FUN = function(x) which(x[,"n_divergent__"] > 0))




stan_model = stan('Gryphon_MMf2.stan',chains=0)

A = 2 * kinship2::kinship(Ped[,1],Ped[,2],Ped[,3])

# a = matrix(rnorm(nrow(A)*10),nc=10);a=a %*% t(a) + diag(10,ncol(A))
# colnames(a) = rownames(a) = colnames(A)
# A = diag(1/sqrt(diag(a))) %*% a %*% diag(1/sqrt(diag(a)))
# colnames(A) = rownames(A) = colnames(a)

i = which(!is.na(Data$BWT))#[1:300]
A_index = match(Data$animal[i],colnames(A))
A_sub = A[A_index,A_index]
X = model.matrix(~SEX,Data[i,])
Z1 = diag(1,ncol(A_sub))
sZAZ = svd(Z1 %*% A_sub %*% t(Z1))
Z2 = model.matrix(~0+BYEAR,Data[i,])
Z3 = model.matrix(~0+MOTHER,Data[i,])
stan_data = list(
    N = length(i),
    Y = Data$BWT[i],
    p_X = ncol(X),
    X = X,
    U = sZAZ$u,
    D = sZAZ$d,
    rZ2 = ncol(Z2),
    rZ3 = ncol(Z3),
    Z2 = Z2,
    Z3 = Z3
    )


Nchains <- 4
Niter <- 200
Nwarm <- 50

stan_fit = sampling(  object = get_stanmodel(stan_model), 
    data =stan_data, 
    iter=Niter,
    warmup = Nwarm,
    init='random',
    chains = Nchains, verbose=TRUE,
    refresh = 10,
    pars=c('beta','sigma2','sigma_Z12','sigma_Z22','sigma_Z32') #
    )
plot(stan_fit)
print(stan_fit)
pairs(stan_fit,pars=c('h2','sigma_p'))
e = do.call(cbind,extract(stan_fit,pars=c('sigma','sigma_Z1','h2')))
h2 = e[,2]^2/(e[,2]^2+e[,1]^2)
plot(h2,e[,3])
abline(0,1)

p.var <- var(Data$BWT, na.rm = TRUE)
prior1.4 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1,
    n = 0.002), G3 = list(V = 1, n = 0.002)), R = list(V = 1,
    n = 0.002))
model1.4 <- MCMCglmm(BWT ~ SEX, random = ~animal + BYEAR + MOTHER,
    pedigree = Ped, data = Data, nitt = 65000, thin = 50, burnin = 15000,
    prior = prior1.4, verbose = T)
summary(model1.4)
posterior.mode(model1.4$VCV)


