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

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("mainFunctions.R")

number.cores <- 16

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm"

## --------------------------------------

xmin <- 11 ; xmax <- 22 ; ymin <- -38 ; ymax <- -21 

resolution <- 0.01

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

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

Meridional.Current

u <- raster(Meridional.Current[2])
v <- raster(Zonal.Current[2])

v <- crop(v,extent(region.as.raster))
u <- crop(u,extent(region.as.raster))

u <- resample(u,region.as.raster,method="bilinear")
v <- resample(v,region.as.raster,method="bilinear")

direction <- 180 + 180 * atan2(v,u) / pi
speed <- sqrt(v^2 + u^2)

# polar = CRS("+init=epsg:3031") # North Pole polar = CRS("+init=epsg:3995")
# direction <- raster::projectRaster(direction, crs = polar, method = "ngb")
# speed <- raster::projectRaster(speed, crs = polar, method = "ngb")

direction[is.na(direction)] <- 0
speed[is.na(speed)] <- 0

currents <- stack(direction,speed)
names(currents) <- c("wind.direction","wind.speed")
plot(currents)

## ------------------

Conductance.min.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS)
Conductance.min.active <- geoCorrection(Conductance.min.active, type="c", multpl=FALSE)

Conductance.mean.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS) 
Conductance.mean.active <- geoCorrection(Conductance.mean.active, type="c", multpl=FALSE)

Conductance.max.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS)
Conductance.max.active <- geoCorrection(Conductance.max.active, type="c", multpl=FALSE)

# Conductance.min.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 
# Conductance.mean.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 
# Conductance.max.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 
# 
# Conductance.min.passive <- geoCorrection(Conductance.min.passive, type="c", multpl=FALSE)
# Conductance.mean.passive <- geoCorrection(Conductance.mean.passive, type="c", multpl=FALSE)
# Conductance.max.passive <- geoCorrection(Conductance.max.passive, type="c", multpl=FALSE)

## ---------------------------------------------
## Distance

cost.surface <- u
cost.surface[!is.na(cost.surface)] <- 1
cost.surface[is.na(cost.surface)] <- 0

Conductance.distance <- transition(cost.surface, mean, directions=8)
Conductance.distance <- geoCorrection(Conductance.distance, type="c", multpl=FALSE)

## ---------------------------------------------
# Thermal

thermal <- raster("/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm/Ocean.temperature.Surface.Var.Lt.Max.tif")
thermal <- crop(thermal,extent(region.as.raster))
thermal <- resample(thermal,region.as.raster,method="bilinear")

speed <- ( 1 / ( thermal / max(getValues(thermal),na.rm=TRUE)) - 1)
speed[is.na(speed)] <- 0

direction <- terrain(thermal, opt=c('flowdir'), unit='degrees')
direction[is.na(direction)] <- 0

plot(speed)
plot(direction)

thermal <- stack(direction,speed)
names(thermal) <- c("wind.direction","wind.speed")

Conductance.thermal <- flow.dispersion(thermal, type ="active", fun=cost.FMGS)
Conductance.thermal <- geoCorrection(Conductance.thermal, type="c", multpl=FALSE)

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

# Isolation by distance (genetic data) vs. Conductance model

project.folder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Drivers of population genetic structurue in Laminaria pallida/Transport Simulation in South Africa/"

genetic.coords.file <- paste0(project.folder,"/Data/Differentiation/Coords.csv")
genetic.diff.file <- paste0(project.folder,"/Data/Differentiation/FST.csv")

genetic.coords <- read.csv(genetic.coords.file,sep=";")[,2:3]
differentiation <- read.csv(genetic.diff.file,sep=";",header=T)[,-1]

differentiation[upper.tri(differentiation)] <- t(differentiation)[upper.tri(differentiation)]

nrow(genetic.coords)
nrow(differentiation)

plot(cost.surface,col=c("#A0CCF2","#737373"),box=FALSE,legend=FALSE)
points(genetic.coords)

## ---------------------------------------------------------

genetic.coords <- relocate.coordinates.na(genetic.coords,cost.surface,maximum.distance=50)
nrow(genetic.coords)

sites <- cellFromXY(cost.surface,genetic.coords)
unique.sites <- unique(sites)

if( length(sites) == length(unique.sites) ) { differentiation <- genetic.diff }

if( length(sites) != length(unique.sites)) {
  
  differentiation.i <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  
  for( i in 1:length(unique.sites)) {
    for( j in 1:length(unique.sites)) {
      
      pairs.i <- which( sites == unique.sites[i]  )
      pairs.j <- which( sites == unique.sites[j]  )
      
      differentiation.i[i,j] <- mean(unlist(differentiation[pairs.i,pairs.j]),na.rm=T)
      
    }  }

  genetic.coords <- genetic.coords[-which(duplicated(sites)),]
  differentiation <- differentiation.i
  
}

diag(differentiation) <- NA

plot(cost.surface,col=c("white","#75C2D3"))
points(genetic.coords,col="#565656",pch=19)
lines( shortestPath(Conductance.distance, as.matrix(genetic.coords[1,]) , as.matrix(genetic.coords[13,]) , output="SpatialLines") )

## ---------------------------------------------------------

# coordinates(genetic.coords) <- ~Lon+Lat
# crs(genetic.coords) <- CRS("+proj=longlat +datum=WGS84")
# summary(genetic.coords)
# 
# genetic.coords <- spTransform(genetic.coords, crs(polar))
# summary(genetic.coords)
# 
# plot(direction)
# points(genetic.coords)
# 
# genetic.coords <- as.data.frame(genetic.coords)

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

  res <- as.numeric(costDistance(Conductance.distance, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))
  if(res == 0) { res <- spDistsN1( as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) , longlat=TRUE) * 1000 }
  return( res )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.currentsDistance.min.active <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.min.active, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)


## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.currentsDistance.max.active <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.max.active, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.currentsDistance.mean.active <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.mean.active, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)


## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.thermalDistance <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.thermal, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## -------------------------------------------

results <- data.frame(comb.all,
                      thermalDistance=results.thermalDistance,
                      geoDistance=results.geoDistance,
                      currentsDistance.min.active=results.currentsDistance.min.active,
                      currentsDistance.mean.active=results.currentsDistance.mean.active,
                      currentsDistance.max.active=results.currentsDistance.max.active,
                      diff=results.diff)

results[results == Inf] <- NA
head(results)

## ----------------

norm <- results[,1:2]
norm <- norm[!duplicated(t(apply(norm, 1, sort))),]

results.final <- data.frame()

for( i in 1:nrow(norm)) {
  
  t.1 <- which( results[,1] == norm[i,1] & results[,2] == norm[i,2] )
  t.2 <- which( results[,1] == norm[i,2] & results[,2] == norm[i,1] )
  
  results.final <- rbind(results.final,
                              data.frame(from=norm[i,1],
                                         to=norm[i,2],
                                         Differantiation = as.numeric(as.character(results$diff[t.1])),
                                         Distance = results$geoDistance[t.1],
                                         
                                         Connectivity.min = min(c(results$currentsDistance.min.active[t.1],results$currentsDistance.mean.active[t.2])) ,
                                         Connectivity.mean = mean(c(results$currentsDistance.mean.active[t.1],results$currentsDistance.mean.active[t.2])) ,
                                         Connectivity.max = max(c(results$currentsDistance.max.active[t.1],results$currentsDistance.mean.active[t.2])) ,
                                         thermalDistance = mean(c(results$thermalDistance[t.1],results$thermalDistance[t.2]),na.rm=T) ,
                                         
                                         stringsAsFactors = FALSE ) )
  
}

## ---------------

results.final[,3] <- results.final [,3] / (1 - results.final [,3])
results.final[results.final[,3] < 0,3] <- 0

results.final[,5] <- sqrt(results.final [,5])
results.final[,6] <- sqrt(results.final[,6] )
results.final[,7] <- sqrt(results.final[,7] )

fit.ibd <- lm(Differantiation ~ Distance, data=results.final , na.action = na.omit)
fit.min.active <- lm(Differantiation ~ Connectivity.min, data=results.final , na.action = na.omit)
fit.mean.active <- lm(Differantiation ~ Connectivity.mean, data=results.final , na.action = na.omit)
fit.max.active <- lm(Differantiation ~ Connectivity.max, data=results.final , na.action = na.omit)
fit.thermal <- lm(Differantiation ~ thermalDistance, data=results.final , na.action = na.omit)

data.frame(aic.ibd= AIC(fit.ibd),
           r2.fit.ibd = summary(fit.ibd)$adj.r.squared,
           fit.min.active= AIC(fit.min.active),
           r2.fit.min.active = summary(fit.min.active)$adj.r.squared,
           fit.mean.active= AIC(fit.mean.active),
           r2.fit.mean.active = summary(fit.mean.active)$adj.r.squared,
           fit.max.active= AIC(fit.max.active),
           r2.fit.max.active = summary(fit.max.active)$adj.r.squared,
           cor.thermal= AIC(fit.thermal),
           r2.fit.thermal = summary(fit.thermal)$adj.r.squared
           
           ) 

plot(results.final$Distance/1000,results.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="Marine distance (km)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results.final$Distance),to=max(results.final$Distance),length.out=100),predict(fit.ibd, data.frame(Distance=seq(from=min(results.final$Distance),to=max(results.final$Distance),length.out=100) )),lty=2,col="#902828")

plot(results.final$Connectivity.mean,results.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="Ocean Currents (cost to disperse)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results.final$Connectivity.mean,na.rm=T),to=max(results.final$Connectivity.mean,na.rm=T),length.out=100),predict(fit.mean.active, data.frame(Connectivity.mean=seq(from=min(results.final$Connectivity.mean,na.rm=T),to=max(results.final$Connectivity.mean,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(results.final$thermalDistance,results.final$Differantiation,lty=1,col="#5E5E5E",ylab="",xlab="Thermal stress (cost to disperse)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results.final$thermalDistance,na.rm=T),to=max(results.final$thermalDistance,na.rm=T),length.out=100),predict(fit.thermal, data.frame(thermalDistance=seq(from=min(results.final$thermalDistance,na.rm=T),to=max(results.final$thermalDistance,na.rm=T),length.out=100) )),lty=2,col="#902828")

## ---------------

fit.mix <- lm(Differantiation ~ thermalDistance+Connectivity.min+Connectivity.mean+Connectivity.max+Distance, data=results.final)
step(fit.mix, direction = c("both"))

fit.final <- lm(Differantiation ~ Connectivity.min+thermalDistance, data=results.final)
AIC(fit.final)
summary(fit.final)
cor(results.final$Differantiation ,predict(fit.final))

## ---------------

par(mfrow=c(1,1),mar = c(6, 6, 2, 2))

plot(results.final$Distance/1000,results.final$Differantiation,lty=1,pch=16,col="#B3B3B3",ylab="",xlab="Marine distance (km)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results.final$Distance),to=max(results.final$Distance),length.out=100),predict(fit.ibd, data.frame(Distance=seq(from=min(results.final$Distance)/1000,to=max(results.final$Distance)/1000,length.out=100) )),lty=2,col="#902828")

par(mfrow=c(1,3),mar = c(5, 5, 4, 4))

plot(results.final$Connectivity.min,results.final$Differantiation,lty=1,pch=19,col="#B3B3B3",ylab="",xlab="Ocean Currents (cost to disperse)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results.final$Connectivity.min,na.rm=T),to=max(results.final$Connectivity.min,na.rm=T),length.out=100),predict(fit.min.active, data.frame(Connectivity.min=seq(from=min(results.final$Connectivity.min,na.rm=T),to=max(results.final$Connectivity.min,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(results.final$thermalDistance,results.final$Differantiation,lty=1,pch=19,col="#B3B3B3",ylab="",xlab="Thermal stress (cost to disperse)",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Genetic differentiation",mgp=c(4,1,0)) 
lines(seq(from=min(results.final$thermalDistance,na.rm=T),to=max(results.final$thermalDistance,na.rm=T),length.out=100),predict(fit.thermal, data.frame(thermalDistance=seq(from=min(results.final$thermalDistance,na.rm=T),to=max(results.final$thermalDistance,na.rm=T),length.out=100) )),lty=2,col="#902828")

plot(results.final$Differantiation,predict(fit.final),lty=1,pch=19,col="#B3B3B3",ylab="",xlab="Observed genetic differentiation",axes=FALSE)
axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()
title(ylab="Predicted genetic differentiation",mgp=c(4,1,0)) 
fit.final.line <- lm(Pred ~ Differantiation, data=data.frame(Differantiation=results.final$Differantiation,Pred=predict(fit.final)), na.action = na.omit)
lines(seq(min(results.final$Differantiation),max(results.final$Differantiation),length.out = 100),predict(fit.final.line, data.frame(Differantiation=seq(min(results.final$Differantiation),max(results.final$Differantiation),length.out = 100))) ,lty=2,col="#902828")

## ---------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------