# ----------------------------------------------- #
# Compare benthic latent model results to other methods 
# ----------------------------------------------- #
#
# CJ Brown 2017-09-21


rm(list = ls())
library(purrr)
library(devtools)
library(lme4)
library(jagstools)
library(tidyr)
library(ggplot2)
library(stringr)
library(PlotTools)
library(vegan)

load_all('~/Code/BenthicLatent')

data(lv_input)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

num_levels <- 2
savname <- paste0('BLM_numlv', num_levels,'_v3.RData')
load(savname)
#
# Get results from LV model 
#

smc <- summary(mcout3)
ibeta <- grep("beta",dimnames(smc[[2]])[[1]])
inuwq <- grep("nuwq",dimnames(smc[[2]])[[1]])
spp_loadings <- smc[[2]][ibeta,3]
spplwr <- smc[[2]][ibeta,1]
sppupr <- smc[[2]][ibeta,5]

site_loadings <- smc[[2]][inuwq,3]
sitelwr <- smc[[2]][inuwq,1]
siteupr <- smc[[2]][inuwq,5]

#
# input variables 
#
y <- lv_input$y
sqrty <- sqrt(lv_input$y)
flow <- lv_input$flow
invlogity <- ((lv_input$y)/375) + 0.01
invlogity <- log(invlogity/(1-invlogity))
mindists <- (lv_input$mindists * lv_input$sddist) + lv_input$mndist


#
# MDS method 
#

mds <- metaMDS(sqrty, distance = "bray", k = 2, autotransform = FALSE)
plot(mds)
ordisurf(mds, lv_input$mindists, add = TRUE, col = "blue")
ord.fit <- envfit(mds ~ mindists, perm = 1000)
ord.fit
plot(ord.fit)

par(mfrow = c(1,2))
plot(mds$species[,1], spp_loadings, xlab = "MDS habitat ordination", ylab = "Bayesian habitat loadings", 
     main = "MDS1", pch = 16, ylim = c(-5, 3))
arrows(mds$species[,1], spplwr, mds$species[,1], sppupr, length = 0)
plot(mds$species[,2], spp_loadings, xlab = "MDS habitat ordination", ylab = "Bayesian habitat loadings", 
     main = "MDS2", pch = 16, ylim = c(-5, 3))
arrows(mds$species[,2], spplwr, mds$species[,2], sppupr, length = 0)
text(-2.15, -2.62, "Halimeda algae", pos = 1, cex = 0.8, xpd = NA)
text(0.68, -2.78, "Halimeda algae", pos = 2, cex = 0.8)


par(mfrow = c(1,2))
plot(mds$points[,1], site_loadings, xlab = "MDS site ordination", ylab = "Bayesian site weights", 
     main = "MDS1", pch = 16, ylim = c(-2, 3))
arrows(mds$points[,1], sitelwr, mds$points[,1], siteupr, length = 0)
plot(mds$points[,2], site_loadings, xlab = "MDS site ordination", ylab = "Bayesian site weights", 
     main = "MDS2", pch = 16, ylim = c(-2, 3))
arrows(mds$points[,2], sitelwr, mds$points[,2], siteupr, length = 0)


plot(site_loadings, mds$points[,1])
plot(site_loadings, mds$points[,2])
plot(site_loadings, mds$points[,3])

#
# Distance gradient method 
#
s <- as.vector(1 - vegdist(sqrty, method = "bray"))
d <- as.vector(dist(mindists, method = "euclidean"))
m1 <- glm(s ~ d, family = binomial(link = log))
summary(m1)
plot(d, s)


#
# PCA 
#

pca <- prcomp(invlogity, scale = T)
plot(pca)
biplot(pca)
summary(lm(pca$x[,1] ~ mindists + flow))
summary(lm(pca$x[,2] ~ mindists + flow))
summary(lm(pca$x[,3] ~ mindists + flow))


#
# CCA 
#
ccamod <- cca(y ~ mindists + flow)
plot(ccamod)
sccamod <- summary(ccamod)
plot(spp_loadings, sccamod$species[,1])
abline(h = 0)
abline(v = 0)

sccamod$concont #82% on first CCA1
sccamod$cont # only 9.77% on first two CCAs
cor(spp_loadings, sccamod$species[,1]) # most of variation here.

plot(sccamod$species[,1], spp_loadings, xlab = "CCA species loadings", ylab = "Bayesian species loadings", pch = 16, ylim = c(-5, 3))
arrows( sccamod$species[,1], spplwr,  sccamod$species[,1], sppupr, length = 0)
text(0.49, -2.78, "Halimeda algae", pos = 2, cex = 0.8)
text(0.46, -0.04, "Algal assemblage", pos = 2, cex = 0.8, srt = 320)

plot(sccamod$sites[,1], site_loadings, xlab = "CCA site loadings", ylab = "Bayesian site loadings", 
       pch = 16, ylim = c(-3, 3))
arrows( sccamod$sites[,1], sitelwr,  sccamod$sites[,1], siteupr, length = 0)
cor(sccamod$sites[,1], site_loadings) ^ 2

plot(spp_loadings, sccamod$species[,1]) # most of variation here.

plot(sccamod$sites[,1], site_loadings)

par(mfrow = c(1,2))
plot(mindists/1000, lv_input$y[,11], xlab = "Distance to nearest log pond (km)", ylab = "Abundance count", 
     pch = 16, main = "Halimeda algae")
plot(mindists/1000, lv_input$y[,3], xlab = "Distance to nearest log pond (km)", ylab = "Abundance count", 
     pch = 16, main = "Algal assemblage")



library(ggplot2)
dat <- data.frame(d = mindists, HA = lv_input$y[,11], AA = lv_input$y[,3])

ggplot(dat, aes(x = d, y = HA)) + 
    geom_point() + 
    stat_smooth(method = "glm", method.args = list(family = "poisson"))

ggplot(dat, aes(x = d, y = AA)) + 
  geom_point() + 
  stat_smooth(method = "glm", method.args = list(family = "poisson"))
