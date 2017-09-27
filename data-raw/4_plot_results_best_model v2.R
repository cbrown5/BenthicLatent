# ----------------------------------------------- #
# Plot results for best model
# ----------------------------------------------- #
#
# CJ Brown 17 Mar 2017
# v2 cleans up v1

rm(list = ls())
library(purrr)
library(devtools)
library(jagstools)
library(tidyr)
library(ggplot2)
library(stringr)
library(PlotTools)
library(dplyr)
library(forcats)
library(sp)
library(rgdal)
library(raster)


load_all('~/Code/BenthicLatent')

data(lv_input)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

num_levels <- 2
savname <- paste0('BLM_numlv', num_levels,'_v3.RData')
# ---------- #
# HABITAT LOADINGS PLOT
# ---------- #

load(savname)
smc <- summary(mcout3)
ibeta <- grep("beta",dimnames(smc[[2]])[[1]])
isort <- order(smc[[2]][ibeta,3], decreasing = T)
smc[[2]][ibeta,3]
lv_input$habnams
dimnames(lv_input$y)[[2]]

dev.new(width = 10, height = 4)
habitat_loadings_plot(savname, prows = 1, pcols = 3)

# ---------- #
# EFFECTS PLOT
# ---------- #
dev.new(width = 5, height = 4)
plot_effects(savname, 'BLM_numlv3_v3.RData', 'BLM_numlv_2_cumdist_v3.RData')

# ------------ #
# Predictions
# ------------- #

dev.new(width = 13, height = 5)
hab_curve_plot(savname)
