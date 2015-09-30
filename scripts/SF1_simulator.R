# --------- Simulator for sparse factor ----------- #
# Written by -- Dan Runcie
# No hierarchical structure.
# Designed to test different n, p, k, or levels of sparsity

generate_correlated_data = function(
	n,  # number of individuals
	p,	# number of traits per individual
	k,	# number of factors
	non_zero_entries, # either a vector or a scaler giving the proportion of non_zero entries in each factor. If a scalar, expanded into a k-vector
	ideosyncratic_variance_parms = c(2,1), # parameters of inverse gamma distribution for ideosyncratic variances
	return_factors = T,
	return_loadings = T,
	return_ideosyncratic_variances = T
	){

	if(length(non_zero_entries) == 1) non_zero_entries = rep(non_zero_entries,k)
	if(length(non_zero_entries) < k) stop('length of non_zero_entries must equal k or 1')

	Lambda = matrix(0,p,k)
	for(i in 1:k){
		ni = p*non_zero_entries[i] # number of non-zero entries in this column
		Lambda[sample(1:p,ni,replace = FALSE),i] = rnorm(ni) # each non-zero entry is drawn from a unit-normal distribution
	}
	F = matrix(rnorm(n*k),n,k)

	ideosyncratic_variances = 1/rgamma(p,shape = ideosyncratic_variance_parms[1], rate = ideosyncratic_variance_parms[2])	


	Y = F %*% t(Lambda) + matrix(rnorm(n*p,0,rep(sqrt(ideosyncratic_variances),each = n)),n,p)

	if(any(!c(return_factors,return_loadings,return_ideosyncratic_variances))) return(Y)

	out = list()
	out$Y = Y
	if(return_factors) out$F = F
	if(return_loadings) out$Lambda = Lambda
	if(return_ideosyncratic_variances) out$ideosyncratic_variances = ideosyncratic_variances
	return(out)
}