# ----------------------------------------------- #
# Simulated data
# ----------------------------------------------- #
#
# CJ Brown 2017-09-08

#Should I do RMSE too? Will require much longer runs, as have to save more params

rm(list = ls())
library(purrr)
library(devtools)
library(lme4)
library(jagstools)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(stringr)
library(PlotTools)
library(dplyr)
library(bmk)


setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

#
# Vary n
#
load("sim-results/varyn.rda")
nseeds <- sum(datsave[[2]] == datsave[[2]][1])

getbeta  <- function(x, expected = seq(-0.2, 0.2, length.out = 5)){
	ibeta <- grep("beta", row.names(x[[1]]$quantiles))
	x[[1]]$quantiles[ibeta,3] - expected
}
dfbeta <- do.call("rbind", map(datsave[[1]], getbeta))

df <- data.frame(n = datsave[[2]],
		dfbeta) %>%
	gather(param, value, -n) %>%
	group_by(n, param) %>%
	summarize(mg = median(value), lwr = quantile(value, 0.25), upr = quantile(value, 0.75))

df$n <- df$n +  rep(seq(-2, 2, length.out = 5), 4)
dev.new()
ggplot(df, aes(x = n, y = mg, color = param)) +
	geom_point() +
	geom_line() +
	geom_linerange(aes(ymin = lwr, ymax = upr)) +
	geom_hline(aes(yintercept = 0)) +
	geom_hline(aes(yintercept = 0.2)) +
	geom_hline(aes(yintercept = 0.1)) +
	geom_hline(aes(yintercept = -0.1)) +
	geom_hline(aes(yintercept = -0.2)) +
	theme_base()


getgamma <- function(x, thresh = 0.3){
#	helldist <- bmkconverge(x$mcout[[1]], binsize = nrow(x$mcout[[1]])/3)
	#isconv <- sum(helldist>thresh) == 0
	#if(isconv){
		igam <- grep("gamma", row.names(x$smry$quantiles))
		xout <- x$smry$quantiles[igam,3]
			#} else {
		#xout <- NA
		#}
}

df <- data.frame(n = datsave[[2]], gamma = map_dbl(datsave[[1]], getgamma)) #%>%
#	dplyr::group_by(n) %>%
#	summarize(mg = mean(gamma, na.rm = T),nruns = sum(!is.na(gamma)), lwr = quantile(gamma, 0.25, na.rm = T), upr = quantile(gamma, 0.75, na.rm = T))

dev.new()
ggplot(df, aes(x = n, y = gamma)) +
	geom_point()
#+ geom_line() +
	#geom_linerange(aes(ymin = lwr, ymax = upr)) +
	#theme_base()

par(mfrow = c(9,2), mar = c(2,2,1,1))
plot(datsave[[1]][[9]]$mcout, auto.layout = F)
datsave[[1]][[116]]$smry$quantiles
#bmkconverge(datsave[[1]][[116]]$mcout[[1]], binsize = nrow(datsave[[1]][[116]]$mcout[[1]])/2)
#bmkconverge(datsave[[1]][[120]]$mcout[[1]], binsize = nrow(datsave[[1]][[120]]$mcout[[1]])/2)


ggplot(df, aes(x = n, y = mb)) +
	geom_point() +
	geom_linerange(aes(ymin = lwrb, ymax = uprb)) +
	geom_hline(aes(yintercept = 0.2)) +
	theme_base()


getgammaCIs <- function(x, kappa1 = 2){
	igam <- grep("gamma", row.names(x[[1]]$quantiles))
	(x[[1]]$quantiles[igam,1] < kappa1) & (x[[1]]$quantiles[igam,5] > kappa1)
}


df <- data.frame(n = datsave[[2]], gamma = map_dbl(datsave[[1]], getgammaCIs)) %>%
	group_by(n) %>%
	summarize(mg = sum(gamma)/nseeds)

ggplot(df, aes(x = n, y = mg)) +
	geom_point() + geom_line() +
	theme_base() +
	ylim(0,1)



#
# varynl, n
#
load("sim-results/varyn-varynlvfitted.rda")
igam <- grep("gamma", row.names(datsave[[1]][[1]][[1]]$quantiles))
getgamma <- function(x){
	c(x[[1]]$quantiles[igam,5], x[[2]]$quantiles[igam,5], x[[3]]$quantiles[igam,5])
}

nseeds <- sum(datsave[[2]] == datsave[[2]][1])
df <- data.frame(n = datsave[[2]]) %>%
	cbind(do.call("rbind", map(datsave[[1]], getgamma))) %>%
	gather(nlvs, gamma, - n) %>%
	group_by(nlvs, n) %>% summarize(mg = mean(gamma), lwr = quantile(gamma, 0.25), upr = quantile(gamma, 0.75))

df$n <- df$n + rep(c(-1, 0, 1), each = 4)
ggplot(df, aes(x = n, y = mg, color = nlvs)) +
	geom_point() + geom_line() +
	geom_linerange(aes(ymin = lwr, ymax = upr)) +
	geom_hline(aes(yintercept = 5)) +
	theme_base()


getgammaCIs <- function(x, kappa1 = 2){
	igam <- grep("gamma", row.names(x[[1]]$quantiles))
	c((x[[1]]$quantiles[igam,1] < kappa1) & (x[[1]]$quantiles[igam,5] > kappa1),
		(x[[2]]$quantiles[igam,1] < kappa1) & (x[[2]]$quantiles[igam,5] > kappa1),
		(x[[3]]$quantiles[igam,1] < kappa1) & (x[[3]]$quantiles[igam,5] > kappa1)
		)
}
df <- data.frame(n = datsave[[2]]) %>%
	cbind(do.call("rbind", map(datsave[[1]], getgammaCIs))) %>%
	#gather(nlvs, gamma, - n) %>%
	group_by(nlvs, n) %>% summarize(mg = sum(gamma)/nseeds)

df$n <- df$n + rep(c(-1, 0, 1), each = 4)
ggplot(df, aes(x = n, y = mg, color = nlvs)) +
	geom_point() + geom_line() +
	theme_base()


rm(datsave)


#
# Vary nhabs
#
load("sim-results/vary-nhabs.rda")

getgamma <- function(x){
	igam <- grep("gamma", row.names(x[[1]]$quantiles))
	x[[1]]$quantiles[igam,3]
}
df <- data.frame(nhabs = datsave[[2]], gamma = map_dbl(datsave[[1]], getgamma)) %>%
	group_by(nhabs) %>%
	summarize(mg = median(gamma), lwr = quantile(gamma, 0.25), upr = quantile(gamma, 0.75))

ggplot(df, aes(x = nhabs, y = mg)) +
	geom_point() +
	geom_linerange(aes(ymin = lwr, ymax = upr))

#
# Vary gamma
#
load("sim-results/vary-kappa.rda")
df <- data.frame(nhabs = datsave[[3]], gamma = map_dbl(datsave, getgamma)) %>%
	group_by(nhabs) %>%
	summarize(mg = mean(gamma), lwr = quantile(gamma, 0.25), upr = quantile(gamma, 0.75))

ggplot(df, aes(x = nhabs, y = mg)) +
	geom_point() +
	geom_linerange(aes(ymin = mg - sdg, ymax = mg + sdg))
