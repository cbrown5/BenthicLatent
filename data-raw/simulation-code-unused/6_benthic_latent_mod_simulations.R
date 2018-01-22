# ----------------------------------------------- #
# Simulated data
# ----------------------------------------------- #
#
# CJ Brown 2017-09-08

rm(list = ls())
library(purrr)
library(devtools)
library(lme4)
library(jagstools)
library(tidyr)
library(ggplot2)
library(stringr)
library(PlotTools)

load_all('~/Code/BenthicLatent')

data(lv_input)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

# ------------- #
# GENERATE DATA
# ------------- #

#test shrinkage priors for gamma
#x <- seq(0, 10, length.out = 100)
#plot(x, dexp(x, 1), type = "l")

#Function that simulates random data set then fits models with different numbers of lvs.
# First LV is the constrained, rest are unconstrained.
runsims <- function(nsites = 50, seed =42, nhabs = 5, nlvs = 3, kappa1 = 0.5,
	num_levels_fit = 1:3, beta1min = -1, beta1max = 1, beta1sd = NA,
	habmnabund = 10,
	nburn = 30000, niter = 100000, nthin = 30){

	set.seed(seed)
	nsamp <- round(niter/nthin)
	#Model params
	beta1 <- rnorm(nhabs * nlvs, mean= 0, sd = 1) %>%
		matrix(nrow = nhabs)
	if(is.na(beta1min)){
		beta1[,1] <- rnorm(nhabs, mean = 0.5, sd = beta1sd)
		} else {
		beta1[,1] <- seq(beta1min, beta1max, length.out = nhabs)
	}
	meanabund <- rep(habmnabund, nhabs)

	# Spatial field

	x1 <- seq(-2, 2, length.out = nsites)
	lv1 <- rnorm(nsites, mean = x1 * kappa1, sd = 1)
	#unconstrained lvs
	lvuc <- rnorm(nsites * (nlvs - 1), mean = 0, sd = 1) %>%
		matrix(nrow = nsites)
	lvmat <- cbind(lv1, lvuc)

	# Linear predictors - standardised to give mean abundance
	linmult <- t(beta1 %*% t(lvmat))
	int <- log(meanabund) - apply(linmult, 2, mean)
	linpred <- linmult + matrix(rep(int, nsites), nrow = nsites, byrow = T)

	# Observations
	obs <- map_df(data.frame(linpred), ~rpois(nsites, lambda = exp(.x))) %>% as.matrix()
	obs[is.na(obs)] <- max(obs, na.rm = T)
	#
	# Fit models
	#
	mcsave <- NULL
	for (ilv in num_levels_fit){
		 jags3 <- jags.model('2a_poissonmodel_matrix_latent_simtesting.txt',
		 data = list('n' = nsites,
		 'y' = obs,
		 'p' = nhabs,
		 'd' = x1,
		 num.lv = ilv
		 ),
		 n.chains = 1,
		inits = list(gamma = kappa1, beta = beta1[,1]), quiet = TRUE)
		update(jags3, nburn)
		mcout3 <- coda.samples(jags3,
		variable.names=c("beta","gamma"),
		n.iter=niter, thin = nthin, quiet = TRUE)
		smc <- summary(mcout3)

		#Just saving summary.
		mcsave <- c(mcsave, list(mcout = mcout3, smry = smc, x1 = x1, kappa1 = kappa1, beta1 = beta1[,1]))
		print(nsites)
		rm(jags3)
		rm(mcout3)
	}
	return(mcsave)
}

# ---------------
# Run simulations
# ---------------
#debugonce(runsims)
#7.8 minutes for one run with 1 mill mcmcs. #5 minutes wiht 15 habs and 100 000 smaples
system.time(
	dout <- runsims(num_levels_fit = 2, nhabs = 8, beta1min = NA, beta1sd = 1)
)

dout[[2]]$quantiles

igamma <-grep("gamma",row.names(dout[[2]]$quantiles))
plot(dout$x1 * dout$kappa1, dout$x1 * dout[[2]]$quantiles[igamma, 3])
abline(0,1)

bmk::bmkconverge(dout$mcout[[1]], binsize = nrow(dout$mcout[[1]])/2)
par(mar = rep(1, 4), mfrow = c(9,2))
plot(dout[[1]], auto.layout = FALSE)

igamma <-grep("gamma",row.names(dout[[2]]$quantiles))
plot(dout[[1]][,igamma])
plot(dout[[1]][,5])

#
# Vary sample size
#

#TODO: First try above code with fewer habs but same low value of Kappa. Can I get convergence faster? Then modify code and run this...

nseeds <- 3
reps <- c(20, 40, 80)
nreps <- length(reps)
rep_vec <- rep(reps, each = nseeds)
seeds <- rep(round(seq(10, 1000, length.out = nseeds)), nreps)

#25 minutes for 9 reps with 8 habs
system.time(
	datout <- map2(rep_vec, seeds, ~runsims(nsites = .x, seed = .y, num_levels_fit = 2, nhabs = 8, beta1min = NA, beta1sd = 1))
)
#save dat out
datsave <- c(list(datout), list(rep_vec, seeds))
save(datsave, file = "sim-results/varyn.rda")
rm(datout, datsave)


#
# Vary sample size and fitted levels
#

nseeds <- 30
reps <- c(10, 20, 40, 80)
nreps <- length(reps)
rep_vec <- rep(reps, each = nseeds)
seeds <- rep(round(seq(10, 1000, length.out = nseeds)), nreps)

#9 minutes for 100 samples...
#42 minutes for 3 x 4
system.time(
	datout <- map2(rep_vec, seeds, ~runsims(nsites = .x, seed = .y))
)
#save dat out
datsave <- c(list(datout), list(rep_vec, seeds))
save(datsave, file = "sim-results/varyn-varynlvfitted.rda")
rm(datout, datsave)

#
# Vary nhabs #47  minute for 4 x 3
#
nhabs <- c(5, 10, 15, 20)
nhabs_len <- length(nhabs)
nseeds <- 5
nhabs_vec <- rep(nhabs, each = nseeds)
seeds <- rep(round(seq(10, 1000, length.out = nseeds)), nhabs_len)
system.time(
	datout <- map2(nhabs_vec, seeds, ~runsims(seed = .y, beta1min = NA, beta1sd = 1,
		nhabs = .x, num_levels_fit = 2))
)
datsave <- c(list(datout), list(nhabs_vec, seeds))
save(datsave, file = "sim-results/vary-nhabs.rda")
rm(datout, datsave)

#
# Vary kappa #14 minutes for 3 x 3
#
kappa <- c(1, 5, 10)
nkappa_len <- length(kappa)
nseeds <- 3
kappa_vec <- rep(kappa, each = nseeds)
seeds <- rep(round(seq(10, 1000, length.out = nseeds)), nkappa_len)
system.time(
	datout <- map2(kappa_vec, seeds, ~runsims(seed = .y, kappa = .x, num_levels_fit = 2))
)
datsave <- c(list(datout), list(kappa_vec, seeds))
save(datsave, file = "sim-results/vary-kappa.rda")
rm(datout, datsave)
