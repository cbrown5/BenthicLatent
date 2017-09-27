# ----------------------------------------------- #
# Graph results of cross validation
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
library(ggthemes)
library(stringr)
library(PlotTools)
library(dplyr)


load_all('~/Code/BenthicLatent')
data(lv_input)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

nhabs <- ncol(lv_input$y)
nsites <- nrow(lv_input$y)
fname <- "sim-results/cross-validation_v3"
num_levels <- 1:4
rmsedat <- data.frame(matrix(NA, nrow = length(num_levels), ncol = nhabs))
names(rmsedat) <- dimnames(lv_input$y)[[2]]

obsprop <- lv_input$y/375

for (ilv in num_levels){
	load(paste0(fname, ilv,".rda"))
#	predprop <- jkpreds / matrix(rep(rowSums(jkpreds), nhabs), ncol = nhabs,byrow = F)
	# rmsedat[ilv,] <- sqrt(colSums((predprop - obsprop)^2)/nsites)
	for (ihabs in 1:nhabs){
		rmsedat[ilv,ihabs] <- -1*sum(dpois(lv_input$y[,ihabs], lambda = jkpreds[,ihabs], log = T))
	}

}

rmsedat$num_levels <- num_levels
rmseplot <- gather(rmsedat, habitat, rmse, -num_levels)

ppi <- 300
width <- 10
asp <- 0.6

# Rename levels
datnams <- read.csv("Appendix_S3_Benthic_Variables.csv")
rmseplot <- left_join(rmseplot, datnams, by = c("habitat" = "Name.in.data.file"))


png(file = "~/Code/BenthicLatent/ms/cross-validation.png",
    width = width * ppi, height = width*ppi*asp,
    res = ppi, antialias = 'none')

ggplot(rmseplot, aes(x = num_levels, y = rmse))  +
geom_line() +
facet_wrap(~Group.for.analysis, scales = "free") +
xlab("Number of unconstrained latent variables") +
ylab("Negative log likelihood") +
theme_bw()#+
#theme(strip.text.x = element_text(size = 8))
dev.off()
