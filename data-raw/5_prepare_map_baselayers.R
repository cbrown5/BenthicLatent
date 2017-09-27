# ---------------
# Prepare baselayers for mapping
# ---------------
# CJ Brown 16 Jun 2017
# This version just uses inverse probit transform

rm(list = ls())
library(devtools)
library(jagstools)
library(tidyr)
library(PlotTools)
library(dplyr)
library(sp)
library(rgdal)
library(raster)
library(RColorBrewer)
library(png)
library(maptools)
library(purrr)

# ---------------
# Parameters
# ---------------
#Quantiles for area affected
quants <- c(0.95, 0.75, 0.5, 0.25)

# ---------------
# Load data
# ---------------

load_all('~/Code/BenthicLatent')
data(lv_input)
data(gdist)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

sitesdf <- read.csv('JuvUVCSites_with_ReefTypes_16Jun2016.csv')

num_levels <- 2
savname <- paste0('BLM_numlv', num_levels,'_v3.RData')
load(savname)

reefs <- readOGR("Kia reefs", "Reef_Strata_Kia")
land <- readOGR("LandPoly", "LandPoly")
ponds <- readOGR("Kia_Logging_Ponds", "Kia_Logging_Ponds")
landuse <- raster('GRID_sol_update/sol_update')

SI <- readOGR("SI_admin", "SLB_adm1") #shapefile of SI
# see http://www.diva-gis.org/datadown to download this.
SI <- rgeos::gSimplify(SI, 0.01)

narr <- readPNG("north-arrow-2.png")
arrrat <- dim(narr)[1]/dim(narr)[2]

#New log ponds
newponds <- data.frame(rbind(c(-7.43639, 158.19973), c(-7.42830, 158.20321), c(-7.42558, 158.20111)))
coordinates(newponds) <- ~X2 + X1
proj4string(newponds) <- CRS('+init=epsg:4326')
newpondsutm <- spTransform(newponds, crs(land))

#
# Filter and transform reefs
#
juvhab <- c('Lagoon exposed fringing - reef flat', 'Intra-lagoon patch-reef complex - subtidal reef flat')
reefs2 <- spTransform(reefs, CRS(proj4string(ponds))) %>%
    subset(C_L3Att_L4 %in% juvhab)


#plot(reefs2)
#plot(land, col = "wheat")
#plot(ponds, col = 'red', add = T)
#plot(landuse)

# ---------------
# Distance raster
# ---------------

#
# Code to estimate distances
#

#landsea <- landuse
#landsea[landuse ==7 | landuse ==4] = 1
#landsea[landuse !=7 & landuse !=4] = NA
# plot(landsea, maxpixels = 5000)
#ls2 <- aggregate(landsea, fact = 5)
#ls2 <- landsea
#iloc <- cellFromXY(ls2, ponds)
#ls2[iloc] <- 2
#ls3 <- ls2
#iloc2 <- cellFromXY(ls3, newpondsutm)
#ls3[iloc2] <- 2
#takes about 30 seconds for full res, <1 second for aggregation at factor 5
#system.time(gdist2 <- gridDistance(ls3, origin = 2, omit  = NA))
#devtools::use_data(gdist2, pkg = "~/Code/BenthicLatent")

#
# Here's a distane layer prepared earlier
#
data(gdist2)
plot(gdist)
plot(reefs2, add = T)

dev.new()
plot(gdist2)
plot(reefs2, add = T)

rreef <- rasterize(reefs2, gdist)
rreefdist <- mask(gdist, rreef)

plot(rreefdist)

# ---------------
# Distance layer
# ---------------

#Transform distances to same scale as used in model
rX1 <- (gdist - lv_input$mndist)/lv_input$sddist
xfrommodel <- lv_input$mindists
gmax <- max(xfrommodel)
gmin <- min(xfrommodel)
rX1[rX1 > gmax] <- gmax
rX1[rX1 < gmin] <- gmin
plot(rX1)
plot(reefs2, add = T)

#Distances for new log ponds
rX1new <- (gdist2 - lv_input$mndist)/lv_input$sddist
rX1new[rX1new > gmax] <- gmax
rX1new[rX1new < gmin] <- gmin
plot(rX1new)
plot(reefs2, add = T)

#Check Transform
plot(raster::extract(rX1, cbind(sitesdf$coordx, sitesdf$coordy)), xfrommodel)
abline(0,1)

# -------------------
# Bayes model params
# -------------------
colnams <- dimnames(mcout3[[1]])[[2]]
igam <- stringr::str_detect(colnams, '\\bgamma\\b')
igamf <- stringr::str_detect(colnams, 'gammaf')
mcmat <- as.matrix(mcout3[[1]])
gamma1 <-  mcmat[,igam]
gammaf <- mcmat[,igamf]
gmode <- quantile(gamma1, 0.5)
g025 <- quantile(gamma1, 0.025)
g975 <- quantile(gamma1, 0.975)
gfmode <- quantile(gammaf, 0.5)

head(mcout3[[1]])

# ---------------
# Predict latent variable spatially
# ---------------

#
# Conceptual dev
#
dseq <- seq(gmin, gmax, length.out = 10)
px <- 1 - pnorm(gmode * dseq, mean = 0, sd = 1)
px2 <- 1 - pnorm(g025 * dseq, mean = 0, sd = 1)
plot(2 * dseq ,1 - pnorm(2 * dseq, mean = 0, sd = 1), type = 'l', ylim = c(0,1), lwd = 2)
lines(0.1 * dseq, 1 - pnorm(0.1 * dseq, mean = 0, sd = 1), col = "red", lty = 2, lwd = 2)
points(gmode * dseq, px, col = "darkgreen", pch = 16)
#
# Probit transfrom
#

rprobit <- rprobit025 <- rprobit975 <- raster(rX1)
rprobit[] <- 1 - pnorm(gmode * rX1[], mean = 0, sd = 1)
rprobit025[] <- 1 - pnorm(g025 * rX1[], mean = 0, sd = 1)
rprobit975[] <- 1 - pnorm(g975 * rX1[], mean = 0, sd = 1)

rprobit_new <- rprobit_new025 <- rprobit_new975 <- raster(rX1)
rprobit_new[] <- 1 - pnorm(gmode * rX1new[], mean = 0, sd = 1)
rprobit_new025[] <- 1 - pnorm(g025 * rX1new[], mean = 0, sd = 1)
rprobit_new975[] <- 1 - pnorm(g975 * rX1new[], mean = 0, sd = 1)

plot(rprobit)
plot(rprobit025)
plot(rprobit975)

# ---------------
# Crop latent variable.
# ---------------
r2 <- mask(rprobit, rreef)
r2_025 <- mask(rprobit025, rreef)
r2_975 <- mask(rprobit975, rreef)
r2new <- mask(rprobit_new, rreef)
r2new_025 <- mask(rprobit_new025, rreef)
r2new_975 <- mask(rprobit_new975, rreef)

eregion <- extent(393286.9, 455857.4, 9136901, 9189922)
rcrop <- crop(rprobit, eregion)
rcropnew <- crop(rprobit_new, eregion)

# ---------------
# Estimates of area affected
# ---------------
nreefs <- sum(!is.na(r2[]), na.rm = T)
(area1 <- map_dbl(quants, ~sum(r2[] >= .x, na.rm = T)/nreefs))
(area1_025 <- map_dbl(quants, ~sum(r2_025[] >= .x, na.rm = T)/nreefs))
(area1_975 <- map_dbl(quants, ~sum(r2_975[] >= .x, na.rm = T)/nreefs))
#new reefs
(areanew <- map_dbl(quants, ~sum(r2new[] >= .x, na.rm = T)/nreefs))
(areanew_025 <- map_dbl(quants, ~sum(r2new_025[] >= .x, na.rm = T)/nreefs))
(areanew_975 <- map_dbl(quants, ~sum(r2new_975[] >= .x, na.rm = T)/nreefs))

#additional area affected, hectares:
quants
areaha <- ((prod(res(r2new)) * nreefs * areanew) - (prod(res(r2new)) * nreefs * area1)) / 10000
areaha025 <- ((prod(res(r2new)) * nreefs * areanew_025) - (prod(res(r2new)) * nreefs * area1_025)) / 10000
areaha975 <- ((prod(res(r2new)) * nreefs * areanew_975) - (prod(res(r2new)) * nreefs * area1_975)) / 10000

#Save values to use in manuscript...
datareas <- data.frame(quants = quants, area_base = c(area1, area1_025, area1_975),
    area_new = c(areanew, areanew_025, areanew_975), area_ha = c(areaha, areaha025, areaha975))

save(datareas, file = "area_reef_change.rda")

# ---------------
# Rotations
# ---------------
pquant <- as(rprobit, "SpatialPixelsDataFrame")
surv <- data.frame(x = sitesdf$coordx, y= sitesdf$coordy)
pregion <- as(extent(pquant), "SpatialPolygons")
pregion <- SpatialPolygonsDataFrame(pregion, data = data.frame(ID = 1))
proj4string(pregion) <- proj4string(rprobit)
#Upload this to http://www.globalforestwatch.org to get an estimate of forest area
#writeOGR(pregion, dsn ="pregion", layer = "pregion", driver = "ESRI Shapefile")

pregion2 <- spTransform(pregion, proj4string(SI))
# ---------------
# Plot
# ---------------
ncols <- 20
ppi <- 300
width <- 8
asp <- 1
fcols <- colorRampPalette((brewer.pal(9, 'Reds')))(ncols)
cols <- hexalpha(fcols, 0.5)
landcol <- "grey"
reefcol <- "grey20"
pondcol <- "black"
survcol <- "black"
survcolbg <- hexalpha("white", 0.5)

zat <- seq(0, 1.0, length.out = ncols+1)
zlabels <- signif(zat, 2)
atlabs <- seq(1, ncols+1, by = 4) #label locations, use every nth one

png(file = "~/Code/BenthicLatent/ms/latentquantmap.png",
    width = width * ppi, height = width*ppi*asp,
    res = ppi, antialias = 'none')
# dev.new(width = width, height = width * asp)
par(mar = c(3,3,3,3), fig = c(0,1,0,1))
plot(land, col = landcol, border = landcol,lwd = 0.8, xlim = c(407995, eregion[2]), density = 30)
plot(pquant, add = T, what = "image", col = cols, zlim = c(0,1))
plot(reefs2, add = T, col = NA, border = reefcol, lwd = 0.7)
plot(ponds, add = T, pch = 16, col = "black")
plot(newpondsutm, add = T, pch = 16, col = "grey60")
points(surv, pch = 21, col = survcol, bg = survcolbg)

fields::image.plot(zlim = c(1,6), legend.only = T, col = cols, breaks = zat,
    smallplot = c(0.76, 0.8, 0.67, 0.87),
    axis.args = list(at = zat[atlabs], labels = zlabels[atlabs]),
    legend.args = list(text = "Probability", side = 3, line = 0.3))

box()

text(424900, 9169091+1000, "Barora Faa", srt = 330)
text(435069.1, 9160313.48, "Barofa Ite", srt = 330)

legend(433828, 9184800, legend = c("Log ponds","Recent log ponds", "Survey sites", "Land", "Lagoonal reefs"),
    pch = c(16, 16, 21, NA, NA),
    fill = c("white", NA, NA, landcol, NA),
    density = c(NA, NA, NA, 40, NA),
    border = c("white", NA, NA, landcol, reefcol),
    bty = "n",
    col = c("black", "grey30", "black", NA, NA)
    )

raster::scalebar(5000, xy = c(407547, 9140442), type = "line",
    label = c(NA, "5km", NA))
text(406969.2, 9138517, "UTM Zone 57S. Coastlines derived from landsat",
    cex = 0.8, col = "grey60", pos = 4, offset = 0)

xleftA = 408559.8; ybottomA = 9185733; asizeA <- 0.03;
rasterImage(narr, xleft = xleftA, ybottom = ybottomA, xright = xleftA + (eregion[2] - eregion[1])*asizeA,
    ytop = ybottomA + (eregion[2] - eregion[1])*asizeA*arrrat, interpolate = F)

# Inset plot
par(fig = c(0.17, 0.47, 0.1, 0.3), new = T)
plot(SI, lwd = 0.5, col = "grey80", xlim = c(157.19, 162.96),
    ylim = c(-12.35, -5.45), border = "grey70")
points(159.95, -9.43, pch = 16, cex = 0.5)
box()
plot(pregion2, add = T, lwd = 1, lty = 1)
text(158, -7, "Study region", font = 3, cex = 0.8, pos = 4)
text(158, -6, "SOLOMON ISLANDS", cex = 0.8)
text(159.95, -9.43, "Honiara", cex = 0.8, pos=4)

dev.off()
