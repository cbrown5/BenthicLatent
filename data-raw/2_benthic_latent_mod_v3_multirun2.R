# ----------------------------------------------- #
# Fit multinomial model with latent variable
# ----------------------------------------------- #
#
#
# Run v3 for multiple levels
# CJ Brown 4 Jan 2017
#New hab classification

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

num_levels <- 1:4

# ---------- #
# SET UP DATA
# ---------- #

#
# Distance metrics
#
mindists <- lv_input$mindists

#
# Sub-set data
#
habnams <- dimnames(lv_input$y)[[2]]
# isub <- c(1, 21, 13, 22, 19) #test groups
isub <- order(habnams)
lv_input$y <- lv_input$y[,isub] #ordered alphabetically
lv_input$p <- ncol(lv_input$y)
nhabs <- lv_input$p
n <- nrow(lv_input$y)

mindist <- lv_input$mindist
habnams <- dimnames(lv_input$y)[[2]]
#
# Matrix df
#

dat2 <- data.frame(mindist = lv_input$mindist, lv_input$y, site = 1:nrow(lv_input$y))

dat3 <- tidyr::gather(dat2, hab, count, -mindist, -site)
dat3$npts <- lv_input$npts[1]

habs <- unique(dat3$hab)
nhabs <- length(habs)

#Covariates - scaled min distance
flowind <- lv_input$flow

habnams

# ---------- #
# RUN MODEL
# ---------- #

nburn <- 30000
niter <- 1000000
nthin <- 50
nsamp <- round(niter/nthin)


#
# Min distance version
#
# num_levels <- 2:4
for (ilv in num_levels){

	 jags3 <- jags.model('2a_poissonmodel_matrix_latent_v3.txt',
	 data = list('n' = nrow(lv_input$y),
	 'y' =lv_input$y,
	 'p' =nhabs,
	 'd' = lv_input$mindists,
	 'flow' = flowind,
	 num.lv = ilv
	 ),
	 n.chains = 1)

	#Burn in
	update(jags3, nburn)

	#Extract samples #takes 13 minutes for 400 000 samples
	system.time(
	mcout3 <- coda.samples(jags3,
	variable.names=c('int', 'beta', 'alpha', 'sigma',
	'nu', 'nuwq', 'gamma', 'gammaf'),
	n.iter=niter, thin = nthin)
	)

	savname <- paste0('BLM_numlv',ilv,'_v3.RData')
	save(mcout3, file = savname)

	rm(jags3)
	rm(mcout3)
}


#
# No distance version
#

library(boral)

bmod1 <- boral(lv_input$y, family = 'poisson', num.lv = 2,
mcmc.control = list(n.burnin = nburn,
	n.iteration = niter, n.thin =nthin))

savname <- 'BLM_numlv_2_ord_v3_boral.RData'
	save(bmod1, file = savname)

#
# No distance, my version...
#

 jags3 <- jags.model('2b_poissonmodel_matrix_latent_nofixed_v3.txt',
	 data = list('n' = nrow(lv_input$y),
	 'y' =lv_input$y,
	 'p' =nhabs,
	 num.lv = 2
	 ),
	 n.chains = 1)

	#Burn in
	update(jags3, nburn)

	#Extract samples
	system.time(
	mcout3 <- coda.samples(jags3,
	variable.names=c('int', 'alpha', 'sigma',
	'nu'),
	n.iter=niter, thin = nthin)
	)

	savname <- 'BLM_numlv_2_ord_v3.RData'
	save(mcout3, file = savname)
