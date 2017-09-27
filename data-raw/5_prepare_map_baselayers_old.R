# ---------------
# Prepare baselayers for mapping
# ---------------
# CJ Brown 28 March 2017

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
quants <- c(0.5, 0.4, 0.25, 0.05)

# ---------------
# Load data
# ---------------

load_all('~/Code/BenthicLatent')
data(lv_input)
data(gdist)
setwd('/Users/s2973410/Code/BenthicLatent/data-raw')

dists <- read.csv('distances_pondsrows_to_sitescols.csv') %>%
  dplyr::select(-X) %>% as.matrix()#remove row numbers

sitesdf <- read.csv('JuvUVCSites_with_ReefTypes_16Jun2016.csv')

num_levels <- 2
savname <- paste0('BLM_numlv', num_levels,'_v3.RData')
load(savname)

reefs <- readOGR("Kia reefs", "Reef_Strata_Kia")
land <- readOGR("LandPoly", "LandPoly")
ponds <- readOGR("Kia_Logging_Ponds", "Kia_Logging_Ponds")
landuse <- raster('GRID_sol_update/sol_update')

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
landsea <- landuse
landsea[landuse ==7 | landuse ==4] = 1
landsea[landuse !=7 & landuse !=4] = NA
# plot(landsea, maxpixels = 5000)

#Change this to use a finer resolution.
#ls2 <- aggregate(landsea, fact = 5)
ls2 <- landsea
iloc <- cellFromXY(ls2, ponds)
ls2[iloc] <- 2
ls3 <- ls2
iloc2 <- cellFromXY(ls3, newpondsutm)
ls3[iloc2] <- 2

#takes about 30 seconds for full res, <1 second for aggregation at factor 5
system.time(gdist2 <- gridDistance(ls3, origin = 2, omit  = NA))

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

# ---------------
# Conceptualising how to visualise influence gradient...
# ---------------

#x <- seq(0, 4, length.out = 100) #vector of distances for sites
#sdx <- 1 #sd of latent
#a1 <- 2 #weak gradient = more dispersal
#a2 <- 6 #strong gradient = less dispersal
#Plots of probability that each site is worse than a site with a latent value = 0
#dev.new()
#plot(x, pnorm(0, mean = x*a1, sd = sdx), type = 'l', xlab = "distance") #weak gradient
#lines(x, pnorm(0, mean = x*a2, sd = sdx), col = "red", lwd = 2) #strong gradient
#points(x, pnorm(0, mean = x*a1, sd = sdx*0.33), col = "green", pch = 16) #less variation
# a2/a1 = 3 and points() with sd = 1/3 is equivalent to a2. Therefore the sd on the
# constrained latent is interchangable with the slope on the distance. So I am correct in
# fixing the sd = 1.


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
#5 - closest = col 173
#23 - far  = col 191
#6 - close = col 174
#prob worse than the most affected site
sitesdf$probworse <- rep(NA, nrow(sitesdf))
X1 <- rX1[]
iland <- which(is.finite(X1))
X1vals <- X1[iland]
X1new <- rX1new[]

#Latent variable at most affected site (closest site)
latent_close <- gmin * gmode#mcmat[,igam]
lowtailL <- FALSE
#Alternate version with distance (unaffected) site
#latent_close <- 0 * gmode#mcmat[,igam] #choosing zero as it splits north and south region in terms of distance
#lowtailL <- FALSE
nvals <- nrow(mcmat)
# Calculate probability each pixel is worse than the closest site
getprob <- function(x, gamval, lc){
    pnorm(lc, mean = x * gamval, sd = 1, lower.tail = lowtailL) #set lower tail = TRUE
}
#This version integrates across all estimates - very slow
getprob2 <- function(x, gamvals, lc){
    mean(map_dbl(gamvals, ~pnorm(lc, mean = x * .x, sd = 1)))
}


#Use gmode if don't want to integrate across uncertainty in rate of change.
system.time(
    latentx <- map_dbl(X1, ~getprob(.x, gmode, latent_close))
)
#Do 95% C.I.s on gamma
latentx025 <- map_dbl(X1, ~getprob(.x, g025, latent_close))
latentx975 <- map_dbl(X1, ~getprob(.x, g975, latent_close))

latentx_new <- map_dbl(X1new, ~getprob(.x, gmode, latent_close))
latentx_new025 <- map_dbl(X1new, ~getprob(.x, g025, latent_close))
latentx_new975 <- map_dbl(X1new, ~getprob(.x, g975, latent_close))

#This one is slow - integrates across uncertainty in gamma
# In the end , more reefs will be affected, bc gradient is flatter...
#system.time(
#    latentx2 <- map(X1vals, ~getprob2(.x, mcmat[1:1000,igam], latent_close))
#)


rlatent <- raster(rX1)
rlatent[] <- latentx
rlatent025 <- raster(rX1)
rlatent025[] <- latentx025
rlatent975 <- raster(rX1)
rlatent975[] <- latentx975

#rlatent[iland] <- unlist(latentx2)

rlatentnew <- raster(rX1)
rlatentnew[] <-latentx_new
rlatentnew025 <- raster(rX1)
rlatentnew025[] <-latentx_new025
rlatentnew975 <- raster(rX1)
rlatentnew975[] <-latentx_new975

plot(rlatent, col = brewer.pal(9, "Reds"))
dev.new(); plot(rlatentnew, col = brewer.pal(9, "Reds"))

# ---------------
# Crop latent variable.
# ---------------
r2 <- mask(rlatent, rreef)
r2_025 <- mask(rlatent025, rreef)
r2_975 <- mask(rlatent975, rreef)
r2new <- mask(rlatentnew, rreef)
r2new_025 <- mask(rlatentnew025, rreef)
r2new_975 <- mask(rlatentnew975, rreef)
eregion <- extent(393286.9, 455857.4, 9136901, 9189922)
rcrop <- crop(rlatent, eregion)
rcropnew <- crop(rlatentnew, eregion)

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
pquant <- as(rlatent, "SpatialPixelsDataFrame")
surv <- data.frame(x = sitesdf$coordx, y= sitesdf$coordy)

#doesn't quite work. Need extents to be the same. Also turns pixesl into points.
#rot <- -30
#pquantr <- elide(pquant, rotate = rot)
#landr <- elide(land, rotate = rot)
#reefs2r <- elide(reefs2, rotate = rot)
#pondsr <- elide(ponds, rotate = rot)

#coordinates(surv) <- ~x + y
#proj4string(surv) <- proj4string(land)
#survr <- elide(surv, rotate = rot)


# ---------------
# Plot
# ---------------
ncols <- 8
ppi <- 300
width <- 8
asp <- 1
fcols <- brewer.pal(ncols, 'Reds')
cols <- hexalpha(fcols, 0.5)
landcol <- "grey"
reefcol <- "black"
pondcol <- "black"
survcol <- "black"
survcolbg <- hexalpha("white", 0.5)

zat <- seq(floor(rlatent@data@min*10)/10, 0.5, length.out = ncols+1)
zlabels <- signif(zat, 2)
atlabs <- seq(1, ncols+1, by = 2) #label locations, skip every second one

png(file = "~/Code/BenthicLatent/ms/latentquantmap.png",
    width = width * ppi, height = width*ppi*asp,
    res = ppi, antialias = 'none')

par(mar = c(3,3,3,3))
plot(land, col = landcol, border = landcol,lwd = 0.8, xlim = c(407995, eregion[2]), density = 30)
plot(reefs2, add = T, col = reefcol, border = NA)
plot(pquant, add = T, what = "image", col = cols)
plot(ponds, add = T, pch = 16, col = "black")
plot(newpondsutm, add = T, pch = 16, col = "grey60")
points(surv, pch = 21, col = survcol, bg = survcolbg)

fields::image.plot(zlim = c(1,6), legend.only = T, col = cols, breaks = zat,
    smallplot = c(0.76, 0.8, 0.67, 0.87),
    axis.args = list(at = zat[atlabs], labels = zlabels[atlabs]),
    legend.args = list(text = "Probability", side = 3, line = 0.3))

box()
legend(433828, 9184800, legend = c("Log ponds","Recent log ponds", "Survey sites", "Land", "Lagoonal reefs"),
    pch = c(16, 16, 21, NA, NA),
    fill = c("white", NA, NA, landcol, reefcol),
    density = c(NA, NA, NA, 40, NA),
    border = c("white", NA, NA, landcol, "black"),
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

dev.off()
