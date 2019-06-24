## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------
##
## Ocean Currents
## Muffins 'n' Code
## https://github.com/jorgeassis
##
##
## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------

source("mainFunctions.R")

number.cores <- 4

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm"
directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Past LGM/Surface/LongTerm"

## --------------------------------------

# xmin <- -140 ; xmax <- -110 ; ymax <- 59 ; ymin <- 25
xmin <- -180 ; xmax <- 180 ; ymax <- 0 ; ymin <- -90

resolution <- 0.1

region.as.table <- matrix( NA ,nrow= ((ymax-ymin)/resolution) ,ncol= ((xmax-xmin)/resolution) )
region.as.raster <- raster(region.as.table)
extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)
crs(region.as.raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

Meridional.Current <- list.files(directory,pattern="Sea.water.X.velocity",full.names = TRUE)
Meridional.Current <- Meridional.Current[which(!grepl("Range", Meridional.Current))]
Meridional.Current <- Meridional.Current[which(!grepl("Summary", Meridional.Current))]

Zonal.Current <- list.files(directory,pattern="Sea.water.Y.velocity",full.names = TRUE)
Zonal.Current <- Zonal.Current[which(!grepl("Range", Zonal.Current))]
Zonal.Current <- Zonal.Current[which(!grepl("Summary", Zonal.Current))]

u <- raster(Meridional.Current[2])
v <- raster(Zonal.Current[2])

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

w <- brick(u, v)
projection(w) <- CRS("+init=epsg:4326")

w2hd.1 <- raster::resample(w[[1]],region.as.raster, method="bilinear")
w2hd.2 <- raster::resample(w[[2]],region.as.raster, method="bilinear")

w2hd <- brick(w2hd.1, w2hd.2)

plot(w2hd)
vectorplot(w2hd * 10, isField = "dXY", region = FALSE, margin = FALSE, narrows = 10000)

slope <- log(sqrt(w2hd[[1]]^2 + w2hd[[2]]^2)) + 9.02
aspect <- atan2(w2hd[[1]], w2hd[[2]])
vectorplot(w2hd * 10, isField = "dXY", region = slope, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 20)),  narrows = 10000, at = seq(0,9.33,length.out=20))

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

polar = CRS("+init=epsg:3031") # North Pole polar = CRS("+init=epsg:3995")

## ------------------

Meridional.Current

u <- raster(Meridional.Current[4])
v <- raster(Zonal.Current[4])

v <- crop(v,extent(region.as.raster))
u <- crop(u,extent(region.as.raster))

direction <- 180 + 180 * atan2(v,u) / pi
speed <- sqrt(v^2 + u^2)

direction <- raster::projectRaster(direction, crs = polar, method = "ngb")
speed <- raster::projectRaster(speed, crs = polar, method = "ngb")

direction[is.na(direction)] <- 0
speed[is.na(speed)] <- 0

currents <- stack(direction,speed)
names(currents) <- c("wind.direction","wind.speed")
plot(currents)

## ------------------

Conductance.min.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS) 
Conductance.mean.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS) 
Conductance.max.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS)

Conductance.min.active <- geoCorrection(Conductance.min.active, type="c", multpl=FALSE)
Conductance.mean.active <- geoCorrection(Conductance.mean.active, type="c", multpl=FALSE)
Conductance.max.active <- geoCorrection(Conductance.max.active, type="c", multpl=FALSE)

Conductance.min.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 
Conductance.mean.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 
Conductance.max.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 

Conductance.min.passive <- geoCorrection(Conductance.min.passive, type="c", multpl=FALSE)
Conductance.mean.passive <- geoCorrection(Conductance.mean.passive, type="c", multpl=FALSE)
Conductance.max.passive <- geoCorrection(Conductance.max.passive, type="c", multpl=FALSE)

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

# Isolation by distance (genetic data) vs. Conductance model

genetic.diff.file <- "/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Genetic/Fst Final Genetic Data 6L S.csv"
genetic.coords.file <- "/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Genetic/summary coords gis.csv" 

genetic.diff <- read.csv(genetic.diff.file,sep=";",header=F)
genetic.coords <- read.csv(genetic.coords.file,sep=";")
genetic.coords <- genetic.coords[genetic.coords$Lat < 0,c(6,5)]

nrow(genetic.coords)
nrow(genetic.diff)

## ---------------

genetic.coords <- relocate.coordinates.na(genetic.coords,u,maximum.distance=25)
nrow(genetic.coords)

sites <- cellFromXY(u,genetic.coords)
unique.sites <- unique(sites)

if( length(sites) == length(unique.sites) ) { differentiation <- genetic.diff }

if( length(sites) != length(unique.sites)) {
  
  differentiation <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  
  for( i in 1:length(unique.sites)) {
    for( j in 1:length(unique.sites)) {
      
      pairs.i <- which( sites == unique.sites[i]  )
      pairs.j <- which( sites == unique.sites[j]  )
      differentiation[i,j] <- mean(unlist(genetic.diff[pairs.i,pairs.j]),na.rm=T)
      
    }  }

  genetic.coords <- genetic.coords[-which(duplicated(sites)),]

}

diag(differentiation) <- NA

## ---------------

coordinates(genetic.coords) <- ~Lon+Lat

crs(genetic.coords) <- CRS("+proj=longlat +datum=WGS84")
summary(genetic.coords)

genetic.coords <- spTransform(genetic.coords, crs(polar))
summary(genetic.coords)

plot(direction)
points(genetic.coords)

genetic.coords <- as.data.frame(genetic.coords)

## ---------------
## IBD

cost.surface <- raster::projectRaster(u, crs = polar, method = "ngb")
cost.surface[!is.na(cost.surface)] <- 1
cost.surface[is.na(cost.surface)] <- 0

Conductance.distance <- transition(cost.surface, mean, directions=8)
Conductance.distance <- geoCorrection(Conductance.distance, type="c", multpl=FALSE)

plot(cost.surface,col=c("#A0CCF2","#737373"),box=FALSE,legend=FALSE)
lines( shortestPath(raster_tr_corrected, as.matrix(genetic.coords[1,]) , as.matrix(genetic.coords[20,]) , output="SpatialLines") )
costDistance(Conductance.distance, as.matrix(genetic.coords[1,]) , as.matrix(genetic.coords[20,]) )
costDistance(Conductance.min.active, as.matrix(genetic.coords[1,]) , as.matrix(genetic.coords[20,]) )

## ---------------

comb.all <- expand.grid(Pair.from=1:nrow(genetic.coords),Pair.to=1:nrow(genetic.coords))
comb.all <- comb.all[comb.all$Pair.from != comb.all$Pair.to,]                 

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.diff <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
return( as.numeric( differentiation[comb.all[i,1],comb.all[i,2]] ) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.geoDistance <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {

  return(  as.numeric(costDistance(Conductance.distance, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) )) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.currentsDistance.mean.active <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.mean.active, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## -------------------------------------------

results <- data.frame(comb.all,geoDistance=results.geoDistance,currentsDistance.mean.active=results.currentsDistance.mean.active,diff=results.diff)
results[ results == Inf] <- NA

# save(results,file = "ModelNorth.RData")
# load("ModelNorth.RData")

# save(results,file = "ModelSouth.RData")
# load(file = "ModelSouth.RData")

## ---------------

results[,4] <- sqrt(results [,4])
results[,5] <- results [,5] / (1 - results [,5])

results[,5] <- sqrt(results [,5])
results[,6] <- sqrt(results [,6])
results[results[,7] < 0,7] <- 0
results[,7] <- results [,7] / (1 - results [,7])

results[,7] <- log(results [,7])
results[,8] <- log(results [,8])
results[,9] <- log(results [,9])

cor.ibd <- cor(results$geoDistance,results$diff , use = "complete.obs",method="pearson")
fit.ibd <- lm(diff ~ geoDistance, data=results , na.action = na.omit)

cor.mean.active <- cor(results$currentsDistance.mean.active,results$diff , use = "complete.obs",method="pearson")
fit.mean.active <- lm(diff ~ currentsDistance.mean.active, data=results , na.action = na.omit)

cor.min.active <- cor(results$currentsDistance.min.active,results$diff , use = "complete.obs",method="pearson")
fit.min.active <- lm(diff ~ currentsDistance.min.active, data=results)
cor.max.active <- cor(results$currentsDistance.max.active,results$diff , use = "complete.obs",method="pearson")
fit.max.active <- lm(diff ~ currentsDistance.max.active, data=results)

cor.min.passive <- cor(results$currentsDistance.min.passive,results$diff , use = "complete.obs",method="pearson")
fit.min.passive <- lm(diff ~ currentsDistance.min.passive, data=results , na.action = na.omit)
cor.mean.passive <- cor(results$currentsDistance.mean.passive,results$diff , use = "complete.obs",method="pearson")
fit.mean.passive <- lm(diff ~ currentsDistance.mean.passive, data=results)
cor.max.passive <- cor(results$currentsDistance.max.passive,results$diff , use = "complete.obs",method="pearson")
fit.max.passive <- lm(diff ~ currentsDistance.max.passive, data=results , na.action = na.omit)

data.frame(aic.ibd= AIC(fit.ibd),
           cor.ibd = cor.ibd,
           fit.mean.active= AIC(fit.mean.active),
           cor.mean.active = cor.mean.active) 

data.frame(aic.ibd= AIC(fit.ibd),
           cor.ibd = cor.ibd,
           fit.min.active= AIC(fit.min.active),
           cor.min.active = cor.min.active,
           fit.mean.active= AIC(fit.mean.active),
           cor.mean.active = cor.mean.active,
           fit.max.active= AIC(fit.max.active),
           cor.max.active = cor.max.active,
           fit.min.passive= AIC(fit.min.passive),
           cor.min.passive= cor.min.passive,
           fit.mean.passive= AIC(fit.mean.passive),
           cor.mean.passive = cor.mean.passive,
           fit.max.passive= AIC(fit.max.passive),
           cor.max.passive = cor.max.passive) 

plot(results$geoDistance,results$diff,lty=1,col="#5E5E5E",ylab="",xlab="Marine distance (km)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results$geoDistance),to=max(results$geoDistance),length.out=100),predict(fit.ibd, data.frame(geoDistance=seq(from=min(results$geoDistance),to=max(results$geoDistance),length.out=100) )),lty=2,col="#902828")

plot(results$currentsDistance.mean.active,results$diff,lty=1,col="#5E5E5E",ylab="",xlab="-log(Min. probability)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results$currentsDistance.mean.active,na.rm=T),to=max(results$currentsDistance.mean.active,na.rm=T),length.out=100),predict(fit.mean.active, data.frame(currentsDistance.mean.active=seq(from=min(results$currentsDistance.mean.active,na.rm=T),to=max(results$currentsDistance.mean.active,na.rm=T),length.out=100) )),lty=2,col="#902828")

## ---------------

r2.ibd = summary(fit.ibd)$adj.r.squared
p.ibd=summary(fit.ibd)$coefficients
p.ibd= ifelse( nrow(p.ibd) == 2 , p.ibd[2,4] , NA)

r2.mean = summary(fit.mean.active)$adj.r.squared
p.mean=summary(fit.mean.active)$coefficients
p.mean= ifelse( nrow(fit.mean.active) == 2 , p.mean[2,4] , NA)

data.frame(aic.ibd= AIC(fit.ibd),
           r2.ibd = r2.ibd,
           cor.ibd = cor.ibd,
           p.ibd=p.ibd,
           aic.mean= AIC(fit.mean),
           r2.mean = r2.mean,
           cor.mean = cor.mean,
           p.mean=p.mean) 

## ---------------
## ---------------
