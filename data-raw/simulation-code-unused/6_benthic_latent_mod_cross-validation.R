# ----------------------------------------------- #
# Fit multinomial model with latent variable
# ----------------------------------------------- #
#
#
# Run v3 for multiple levels
# CJ Brown 2017-09-15
 # leave one out cross validation

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
nburn <- 5000
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
getobs_pred_hyper <- function(isamp, ilv, nhabs, params, mindists, flow, y){
	abund <- exp((params$beta[isamp,] * (params$gammaf[isamp] * flow + mindists * params$gamma[isamp])) +
			params$int[isamp,])
	yobs <- lv_input$y[i,]
	ypred <- abund
	dat <- data.frame(obs = yobs, pred = ypred, isamp = isamp)
	dat$habs <- attr(lv_input$y,"dimnames")[[2]]
	return(dat)
	}

#
# Min distance version
#

jksave <- NULL
#ilv <- 4
 for (ilv in 1:2){
 #for (ilv in 4){
	savname <- paste0('BLM_numlv',ilv,'_v3.RData')
	load(savname)

	smc <- summary(mcout3)
	rnames <- row.names(smc$quantiles)
	igamma <- grep("gamma", rnames)
	isigma <- grep("sigma", rnames)
	inuwq <- grep("nuwq", rnames)
	inu <- grep("nu\\[", rnames)
	ibeta <- grep("beta", rnames)
	sq <- smc$quantiles

	initslist <- list(gamma = sq[igamma[1],3],
		gammaf = sq[igamma[2],3],
		beta = sq[ibeta,1], sigma = sq[isigma,3])

	jkpreds <- matrix(NA, nrow = n, ncol = nhabs)

	for (i in 1:n){
		print(i)
		initemp <- initslist
		initemp$nuwq <- sq[inuwq,3][-i]
		initemp$nu <- matrix(sq[inu,3], ncol = ilv)
		initemp$nu <- matrix(initemp$nu[-i,], ncol = ilv)

		jags_new <- jags.model('2a_poissonmodel_matrix_latent_v3.txt',
		   	 data = list('n' = nrow(lv_input$y)-1,
		   	 'y' =lv_input$y[-i,],
		   	 'p' =nhabs,
		   	 'd' = lv_input$mindists[-i],
		   	 'flow' = flowind[-i],
		   	 num.lv = ilv
		   	 ),
		   	 n.chains = 1,
			inits = initemp, quiet = TRUE)
		update(jags_new, nburn)
		mcout_new <- coda.samples(jags_new,
			variable.names=c('int', 'beta', 'gamma', 'gammaf'),
			n.iter=5000, thin = 5)

		nsamp <- nrow(mcout_new[[1]])
		params <- extract_params(mcout_new)
		preds <- purrr::map(1:nsamp, ~getobs_pred_hyper(.x, ilv, nhabs, params,	lv_input$mindists[i], flowind[i], lv_input$y[i,]))
		datpreds <- do.call("rbind", preds)
		jkpreds[i,] <- tapply(datpreds$pred, datpreds$habs, mean)
	}
	save(jkpreds, file = paste0('sim-results/cross-validation_v3',ilv,'.rda'))
	jksave <- c(jksave, list(jkpreds))
	rm(jkpreds)
}
# savname <- paste0('sim-results/cross-validation_v3.rda')
# save(jksave, file = savname)

# habcov <- lv_input$y / 375

#for(i in 1:nhabs){
	#print(cor(jksave[[1]][,i], habcov[,i]))
# }#
