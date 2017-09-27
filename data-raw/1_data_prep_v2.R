# ----------------------------------------------- #
# Prep data for latent variable model
# ----------------------------------------------- #
#
# CJ Brown 13 Jun 2017
#
# Aggregate to avoid groups with many zeros

rm(list = ls())
library(dplyr)
library(tidyr)
library(readr)
library(rgdal)
library(raster)

library(devtools)
load_all('~/Code/BenthicLatent/')

setwd('/Users/s2973410/Code/BenthicLatent/data-raw')
#
# Data
#

dat <- read_csv('JuvUVCSites_with_ReefTypes_16Jun2016.csv')

head(dat)
nsites <- nrow(dat)

bendat <- read_csv('BenthicCoverSurveys.csv')

nhabs <- length(unique(bendat$code))

# Site 46, transect 5 is missing three points
bendat[7120,'cover'] <- 14
# filter(bendat, site ==46) %>% data.frame()
bendat$n.pts <- 75
bencodes <- read_csv('benthic_codes.csv') %>% distinct()
write.csv(bencodes, "Benthic_Variables.csv", row.names = F)
#
# Distance to log ponds
#

ponds <- readOGR("Kia_Logging_Ponds", "Kia_Logging_Ponds")
landuse <- raster('GRID_sol_update/sol_update')

landsea <- landuse
landsea[landuse ==7 | landuse ==4] = 1
landsea[landuse !=7 & landuse !=4] = NA
# plot(landsea, maxpixels = 5000)

#Change this to use a finer resolution.
#ls2 <- aggregate(landsea, fact = 5)
ls2 <- landsea
iloc <- cellFromXY(ls2, ponds)
ls2[iloc] <- 2

#takes about 30 seconds for full res, <1 second for aggregation at factor 5
system.time(gdist <- gridDistance(ls2, origin = 2, omit  = NA))

#Buffer points at diam slightly larger than pixel size, to make sure you get a value
dat$mindist <- raster::extract(gdist, cbind(dat$coordx, dat$coordy), buffer = 50, fun = mean, na.rm = T)

#
# Summarise bendat by sites
#

#Check points are correct
bcount <- bendat %>% group_by(site) %>%
summarize(tot = sum(cover))#, totn = sum(n.pts)/nhabs)

data.frame(bcount)


benthic_summary <- bendat %>% group_by(site, code) %>%
summarize(pts = sum(cover)) %>% left_join(bcount) %>%
mutate(cov = pts/tot) %>%
dplyr::select(site, code, pts) %>%
spread(code, pts, fill = 0) %>%
ungroup() %>%
left_join(bcount)

dplyr::select(benthic_summary, -site, - tot) %>% rowSums()


#
# Freq of occurrence, to select variables
#

bentbl <- benthic_summary %>%
  summarize_all(function(x) signif(sum(x>0)/sum(x>=0),2)) %>%
  gather(CODE, freqency) %>%
  left_join(bencodes)

# bentbl %>% filter(freqency<0.6) %>%  data.frame()
 bentbl %>%  data.frame()

#
# Select variables
#
 datsa <- benthic_summary %>%
mutate(branching_dead = ACBD + CBD,
algal_assemblage = DCA + AA,
coral_submassive = ACS + CS,
coral_other = ACD + ACE + CE + CF + CMR + SC,
coral_rare = CHL + CME + CTU + CA,
rock = R + RCK,
unconsolidated = SI + S,
other = SP + ZO
) %>%
dplyr::select(-ACBD, -CBD, -DCA, -AA, - ACS, -CS, -ACD, -ACE, -CE, -CF,-CHL,-CMR, -SC, -CHL, -CME, -CTU, -CA, -R, -RCK, -SI,-S, -SP, -ZO)

dplyr::select(datsa, -site, -tot) %>% rowSums() #should all be the same = 375

#
# Input to model
#

#Centre y - matrix of habitats at sites
y <- dplyr::select(datsa, -site, -tot) %>% as.matrix()
npts <- as.numeric(benthic_summary$tot)

#params for variable sizes
N <- nrow(y)
K <- ncol(y)
J <- length(ponds)

#jags inputs
mindists <- (dat$mindist - mean(dat$mindist))/sd(dat$mindist)
flowind <- as.numeric(factor(dat$flow))-1

lv_input <- list(
  npts = npts,
  N = N,
  K = K,
  J = J,
  p = ncol(y),
  y = y,
  flow = flowind,
  mindists = mindists,
  mndist = mean(dat$mindist),
  sddist = sd(dat$mindist)
  )

isub <- order(habnams)
lv_input$y <- lv_input$y[,isub] #ordered alphabetically
habnams <- dimnames(lv_input$y)[[2]]
lv_input$habnams <- habnams

#Make names for labels
habnams_full <- data.frame(CODE = habnams) %>%
    left_join(bencodes)
habnams_full$CATEGORY[habnams_full$CODE == 'branching_dead'] <- 'Dead branching coral'
habnams_full$CATEGORY[habnams_full$CODE == 'coral_rare'] <- 'Rare corals'
habnams_full$CATEGORY[habnams_full$CODE == 'other'] <- 'Other habitats'
habnams_full$CATEGORY[habnams_full$CODE == 'unconsolidated'] <- 'Soft sediment'
habnams_full$CATEGORY[habnams_full$CODE == 'coral_submassive'] <- 'Submassive coral'
habnams_full$CATEGORY[habnams_full$CODE == 'algal_assemblage'] <- 'Algal assemblage'
habnams_full$CATEGORY[habnams_full$CODE == 'coral_other'] <- 'Other corals'
habnams_full$CATEGORY[habnams_full$CODE == 'rock'] <- 'Rock'

lv_input$habnams_full <- habnams_full

#
# Save data
#

 setwd('../..')

devtools::use_data(lv_input, pkg = 'BenthicLatent', overwrite = T)
devtools::use_data(benthic_summary, pkg = 'BenthicLatent', overwrite = T)
devtools::use_data(bentbl, pkg = 'BenthicLatent', overwrite = T)
devtools::use_data(gdist, pkg = 'BenthicLatent', overwrite = T)

sitesdf <- dat
devtools::use_data(sitesdf, pkg = 'BenthicLatent', overwrite = T)
