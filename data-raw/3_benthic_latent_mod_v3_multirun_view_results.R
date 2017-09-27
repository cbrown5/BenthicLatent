# ----------------------------------------------- #
# View results of model fits
# ----------------------------------------------- #
#
# Run v3 for multiple levels
# CJ Brown 14 Jun 2017
#New hab classification
# This script compares WAIC. It include min distance, cumulativ distance, no distance and boral no distance
# The boral no distance is not directly comparable to my no distance model (2b_poissonmodel_matrix...R) b/c my version had random species intercepts.


rm(list = ls())
library(purrr)
library(devtools)
library(jagstools)
library(tidyr)
library(ggplot2)
library(stringr)
library(PlotTools)
library(bmk)

load_all('~/Code/BenthicLatent')

data(lv_input)
data(sitesdf)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

num_levels <- 1:4
nextra <- 1 #extra draws for marginal LL
dat <- read_csv('JuvUVCSites_with_ReefTypes_16Jun2016.csv')

#
# Functions
#

getlogLik <- function(isamp, ilv, n, nhabs, params){
	abundmat <- mapabund(ilv, n, nhabs, params$nuwq[isamp,],
		 params$nu[isamp,], params$alpha[isamp,], params$beta[isamp,],
		 params$int[isamp,]) %>%
		 unlist() %>% matrix(nrow = n, byrow = T)

	condLogLik <- dpois(lv_input$y, abundmat, log = T)
	apply(condLogLik, 1, sum) #sum across spp, as in Boral
	}

getlogLik_hyper <- function(isamp, ilv, n, nhabs, params){
	abundmat <- mapabund_hyper(ilv, n, nhabs, lv_input$mindists, lv_input$flow, params$alpha[isamp,], params$beta[isamp,],params$int[isamp,],
	params$gamma[isamp], params$gammaf[isamp])  %>%
		 unlist() %>% matrix(nrow = n, byrow = T)

	LogLik <- dpois(lv_input$y, abundmat, log = T)
	apply(LogLik, 1, sum) #sum across spp, as in Boral
	}


rootogram <- function(yobs, ypred, nbrks = 10, ...){
	 brks <- seq(0, ceiling(max(c(ypred, yobs))), length.out = nbrks)
	  xwd <- diff(brks)[1]/2.1
	 xy <- hist(ypred, breaks = brks, plot = F)
	 xy2 <-hist(yobs, breaks = brks, plot = F)
	ylwr <-  sqrt(xy$counts) - sqrt(xy2$counts)
	ylim <- c(min(ylwr), sqrt(max(xy$counts)))

	 plot(xy$mids, sqrt(xy$counts), lty = 3, col = 'red', type = 'b', ylim = ylim,
	 xlab ='Count', ylab = 'sqrt(Frequency)', ...)
	 rect(xy$mids - xwd, sqrt(xy$counts), xy$mids + xwd,ylwr, col = 'grey')
	lines(xy$mids, sqrt(xy$counts), lty = 3, type = 'b', pch = 16, col = 'red')
	abline(h = 0)
	}

	#Calculate variance and bias
varbias <- function(i, abund_mod, abund_obs){
	ca <- cor(abund_mod[,i], abund_obs[,i])
	bias <- as.numeric(coef(lm(abundmat[,i] ~ 0 + lv_input$y[,i]))[1])
	c(ca, bias)
	}



# ---------- #
# SET UP DATA
# ---------- #

#Covariates
x <- lv_input$mindists
flowind <- lv_input$flow
n <- nrow(lv_input$y)
nhabs <- ncol(lv_input$y)
habnams <- dimnames(lv_input$y)[[2]]

meanabund <- as.numeric(apply(lv_input$y, 2, mean))
#
# View results
#

nlvruns <- length(num_levels)
waic <- rep(NA, nlvruns)
lpd <- rep(NA, nlvruns)
pwaic <- rep(NA, nlvruns)

waic_marg <- rep(NA, nlvruns)
lpd_marg <- rep(NA, nlvruns)
pwaic_marg <- rep(NA, nlvruns)

est_varbias <- NULL
hellconv <- NULL
secchi_cor <- NULL

for (iruns in 1:nlvruns){
	ilv <- num_levels[iruns]
	savname <- paste0('BLM_numlv',ilv,'_v3.RData')
	load(file = savname)

	#Check convergence
	#hellconv <- c(hellconv, list(bmkconverge(mcout3[[1]], binsize = nrow(mcout3[[1]])/10)))
	#Values near zero indicate convergence

	lvs <- jagsresults(x = mcout3, params = c('nu'))
	lvswq <- jagsresults(x = mcout3, params = c('nuwq'))

	gamma <- jagsresults(x = mcout3, params = c('gamma'))
	gammaf <- jagsresults(x = mcout3, params = c('gammaf'))


	int <- jagsresults(x = mcout3, params = c('int'))
	alpha <- jagsresults(x = mcout3, params = c('alpha'))
	beta <- jagsresults(x = mcout3, params = c('beta'))

	abundmat <- mapabund(ilv, n, nhabs, lvswq[,1], lvs[,1], alpha[,1], beta[,1], int[,1]) %>%
	unlist() %>% matrix(nrow = n, byrow = T)


	  # for (i in 1:nhabs){
		  # print(cor(lv_input$y[,i], abundmat[,i]))

		 # dev.new()
		 # rootogram( abundmat[,i], lv_input$y[,i], main = habnams[i], nbrks = 8)
		 # # plot(lv_input$y[,i], abundmat[,i], main = habnams[i])
		 # # abline(0,1)
		  # }


	est_varbias <- c(est_varbias,
	list(
	map(1:nhabs, ~varbias(.x, lv_input$y, abundmat)) %>% unlist() %>%
	matrix(nrow = nhabs, byrow = T) %>% data.frame() %>% set_names(c('r2','bias')) %>% cbind(habnams) %>%
	cbind(meanabund)
	)
	)

	params <- extract_params(mcout3)
	nsamp <- nrow(params$alpha)

	print(iruns)

	#
	# Conditional WAIC
	#
	logLikComp <- purrr::map(1:nsamp, ~getlogLik(.x, ilv, n, nhabs, params)) %>%
	purrr::as_vector() %>% matrix(nrow = nsamp, byrow = T)

	lpd_i <- log(apply(exp(logLikComp), 2, mean, na.rm = TRUE))
	pwaic_i <- apply(logLikComp, 2, var, na.rm = TRUE)
	waic[iruns] <- -2*(sum(lpd_i) - sum(pwaic_i))
	lpd[iruns] <- sum(lpd_i)
	pwaic[iruns] <- sum(pwaic_i)

	#
	# marginal waic
	#

	logLikmarg <- purrr::map(1:nsamp, ~getlogLik_hyper(.x, ilv, n, nhabs, params)) %>%
	purrr::as_vector() %>% matrix(nrow = nsamp, byrow = T)

	lpd_marg_i <- log(apply(exp(logLikmarg), 2, mean, na.rm = TRUE))
	pwaic_marg_i <- apply(logLikmarg, 2, var, na.rm = TRUE)
	waic_marg[iruns] <- -2*(sum(lpd_marg_i) - sum(pwaic_marg_i))
	lpd_marg[iruns] <- sum(lpd_marg_i)
	pwaic_marg[iruns] <- sum(pwaic_marg_i)

	summary(mcout3)
	secchi_cor <- c(secchi_cor, list(cor.test(lvswq[,5], sitesdf$secchi)))

}

save(est_varbias, file = 'LV mod estimtes of variance and bias.RData')
save(secchi_cor, file = "Secchi_correlation.rda")
plot(est_varbias[[4]]$meanabund, est_varbias[[4]]$r2)

waic
which.min(waic)
waic_marg
which.min(waic_marg)
# which.min(waic_marg)
pwaic


load('BLM_numlv_2_cumdist_v3.RData')

	lvs <- jagsresults(x = mcout3, params = c('nu'))
	lvswq <- jagsresults(x = mcout3, params = c('nuwq'))

	gamma <- jagsresults(x = mcout3, params = c('gamma'))
	gammaf <- jagsresults(x = mcout3, params = c('gammaf'))


	int <- jagsresults(x = mcout3, params = c('int'))
	alpha <- jagsresults(x = mcout3, params = c('alpha'))
	beta <- jagsresults(x = mcout3, params = c('beta'))

	abundmat <- mapabund(ilv, n, nhabs, lvswq[,1], lvs[,1], alpha[,1], beta[,1], int[,1]) %>%
	unlist() %>% matrix(nrow = n, byrow = T)

	params <- extract_params(mcout3)
	nsamp <- nrow(params$alpha)

	logLikComp <- purrr::map(1:nsamp, ~getlogLik(.x, 2, n, nhabs, params)) %>%
	purrr::as_vector() %>% matrix(nrow = nsamp, byrow = T)

	lpd_i <- log(apply(exp(logLikComp), 2, mean, na.rm = TRUE))
	pwaic_i <- apply(logLikComp, 2, var, na.rm = TRUE)
	 (waic_cumdist <- -2*(sum(lpd_i) - sum(pwaic_i)))


#
# No distance version
#

library(boral)

load('BLM_numlv_2_ord_v3_boral.RData')
bmod1$ics

#boral(lv_input$y, family = 'poisson', num.lv = 2,
#do.fit = F, model.name = "boral_poisson.txt")

#
# No distance version - CB code
#
# This differs from the boral model, b/c I have random intercepts (due to using poisson to model proportions)
#

load('BLM_numlv_2_ord_v3.RData')

lvs <- jagsresults(x = mcout3, params = c('nu'))

int <- jagsresults(x = mcout3, params = c('int'))
alpha <- jagsresults(x = mcout3, params = c('alpha'))

abundmat <- mapabund(ilv, n, nhabs, 0, lvs[,1], alpha[,1], 0, int[,1]) %>%
	unlist() %>% matrix(nrow = n, byrow = T)

params <- extract_params(mcout3)
params$nuwq <- matrix(0, nrow = nrow(mcout3[[1]]), ncol = n)
params$gamma <- matrix(0, nrow = nrow(mcout3[[1]]))
params$gammaf <- matrix(0, nrow = nrow(mcout3[[1]]))
params$beta <- matrix(0, nrow = nrow(mcout3[[1]]), ncol = nhabs)
nsamp <- nrow(params$alpha)

isamp <- 1
abundmat <- mapabund(2, n, nhabs, params$nuwq[isamp,],
			 params$nu[isamp,], params$alpha[isamp,], params$beta[isamp,],
			 params$int[isamp,])

abundmod(49,2, n, nhabs, params$nuwq[isamp,], params$nu[isamp,], params$alpha[isamp,], params$beta[isamp,], params$int[isamp,])

logLikComp <- purrr::map(1:nsamp, ~getlogLik(.x, 2, n, nhabs, params)) %>%
	purrr::as_vector() %>% matrix(nrow = nsamp, byrow = T)

lpd_i_nodist <- log(apply(exp(logLikComp), 2, mean, na.rm = TRUE))
pwaic_i_nodist <- apply(logLikComp, 2, var, na.rm = TRUE)

(waic_nodist <- -2*(sum(lpd_i_nodist) - sum(pwaic_i_nodist)))

logLikmarg <- purrr::map(1:nsamp, ~getlogLik_hyper(.x, 2, n, nhabs, params)) %>%
purrr::as_vector() %>% matrix(nrow = nsamp, byrow = T)

lpd_marg_i <- log(apply(exp(logLikmarg), 2, mean, na.rm = TRUE))
pwaic_marg_i <- apply(logLikmarg, 2, var, na.rm = TRUE)
waic_nodist_marg <- -2*(sum(lpd_marg_i) - sum(pwaic_marg_i))

#
# Save results
#
waic_dat <- list(waic = waic, waic_marg = waic_marg, waic_cumdist = waic_cumdist, waic_nodist = waic_nodist)
save(waic_dat, file = "waic.rda")
