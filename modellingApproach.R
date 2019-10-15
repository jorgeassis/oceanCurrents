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

number.cores <- 8

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers HR/Present/Surface/LongTerm PerSeason"

## --------------------------------------

xmin <- -6.5 ; xmax <- 2.6 ; ymin <- 46.32  ; ymax <- 51.9
resolution <- 0.005

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

u <- raster(Meridional.Current[4])
v <- raster(Zonal.Current[4])

v <- crop(v,extent(region.as.raster))
u <- crop(u,extent(region.as.raster))

u <- resample(u,region.as.raster,method="bilinear")
v <- resample(v,region.as.raster,method="bilinear")

direction <- 180 + 180 * atan2(v,u) / pi
speed <- sqrt(v^2 + u^2)

## ------------------

polar.1 <- CRS("+init=epsg:3031") # North Pole polar = CRS("+init=epsg:3995")
polar.2 <- CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
direction <- raster::projectRaster(direction, crs = polar.2, method = "ngb")
speed <- raster::projectRaster(speed, crs = polar.2, method = "ngb")

## ------------------

direction[is.na(direction)] <- 0
speed[is.na(speed)] <- 0

currents <- stack(direction,speed)
names(currents) <- c("wind.direction","wind.speed")
plot(currents)

## ------------------

Conductance.mean.active <- flow.dispersion(currents, type ="active", fun=cost.FMGS) 
Conductance.mean.active <- geoCorrection(Conductance.mean.active, type="c", multpl=FALSE)

Conductance.mean.passive <- flow.dispersion(currents, type ="passive", fun=cost.FMGS) 
Conductance.mean.passive <- geoCorrection(Conductance.mean.passive, type="c", multpl=FALSE)

## ---------------------------------------------
## Distance

cost.surface <- u
cost.surface[!is.na(cost.surface)] <- 1
cost.surface[is.na(cost.surface)] <- 0

cost.surface <- raster::projectRaster(cost.surface, crs = polar.2, method = "ngb")

Conductance.distance <- transition(cost.surface, mean, directions=8)
Conductance.distance <- geoCorrection(Conductance.distance, type="c", multpl=FALSE)

## ---------------------------------------------
# Thermal

thermal <- raster("/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers HR/Present/Surface/LongTerm/Ocean.temperature.Surface.Var.Lt.Max.tif")
thermal <- crop(thermal,extent(region.as.raster))
thermal <- resample(thermal,region.as.raster,method="bilinear")
thermal <- thermal - min(getValues(thermal),na.rm=TRUE)
thermal <- thermal / max(getValues(thermal),na.rm=TRUE)

direction <- terrain(thermal, opt=c('aspect'), unit='degrees')
speed <- terrain(thermal, opt=c('slope'), unit='degrees')

# direction <- raster::projectRaster(direction, crs = polar.2, method = "ngb")
# speed <- raster::projectRaster(speed, crs = polar.2, method = "ngb")

speed[is.na(speed)] <- 0
direction[is.na(direction)] <- 0

plot(speed)
plot(direction)

thermalFlow <- stack(direction,speed)
names(thermalFlow) <- c("wind.direction","wind.speed")
plot(thermalFlow)

ConductanceThermalA <- flow.dispersion(thermalFlow, type ="active", fun=cost.FMGS)
ConductanceThermalA <- geoCorrection(ConductanceThermalA, type="c", multpl=FALSE)

ConductanceThermalP <- flow.dispersion(thermalFlow, type ="passive", fun=cost.FMGS)
ConductanceThermalP <- geoCorrection(ConductanceThermalP, type="c", multpl=FALSE)

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

# Isolation by distance (genetic data) vs. Conductance model

project.folder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Genetic differentiation of Laminaria digitata across Brittany/Data/"

genetic.coords.file <- paste0(project.folder,"/39popFinal.csv")
genetic.diff.file <- paste0(project.folder,"/fst.csv")
genetic.diff.file2 <- paste0(project.folder,"/nei.csv")
genetic.diff.file3 <- paste0(project.folder,"/migration.csv")

genetic.coords <- read.csv(genetic.coords.file,sep=";",dec=",")[,2:3]
genetic.coords.names <- read.csv(genetic.coords.file,sep=";",dec=",",stringsAsFactors = F)[,1]
differentiation <- read.csv(genetic.diff.file,sep=";",header=F,dec=",",stringsAsFactors = F)
differentiation2 <- read.csv(genetic.diff.file2,sep=";",header=F,dec=",",stringsAsFactors = F)
differentiation3 <- read.csv(genetic.diff.file3,sep=";",header=F,dec=",",stringsAsFactors = F)
differentiation3[is.na(differentiation3)] <- 0

nrow(genetic.coords)
nrow(differentiation)

## ---------------------------------------------------------

#coordinates(genetic.coords) <- ~Lon+Lat
#crs(genetic.coords) <- CRS("+proj=longlat +datum=WGS84")
#summary(genetic.coords)
 
#genetic.coords <- spTransform(genetic.coords, crs(polar.1))
#summary(genetic.coords)
#genetic.coords <- as.data.frame(genetic.coords)

## ---------------------------------------------------------

plot(cost.surface,col=c("#A0CCF2","#737373"),box=FALSE,legend=FALSE)
points(genetic.coords)

## ---------------------------------------------------------

genetic.coords <- relocate.coordinates.na(genetic.coords,cost.surface,maximum.distance=Inf)
nrow(genetic.coords)

sites <- cellFromXY(cost.surface,genetic.coords)
unique.sites <- unique(sites)

if( length(sites) != length(unique.sites) )  {
  
  differentiation.i <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  differentiation.i2 <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  differentiation.i3 <- matrix(NA,ncol=length(unique.sites),nrow=length(unique.sites))
  
  for( i in 1:length(unique.sites)) {
    for( j in 1:length(unique.sites)) {
      
      pairs.i <- which( sites == unique.sites[i]  )
      pairs.j <- which( sites == unique.sites[j]  )
      
      differentiation.i[i,j] <- mean(unlist(differentiation[pairs.i,pairs.j]),na.rm=T)
      differentiation.i2[i,j] <- mean(unlist(differentiation2[pairs.i,pairs.j]),na.rm=T)
      differentiation.i3[i,j] <- mean(unlist(differentiation3[pairs.i,pairs.j]),na.rm=T)
      
    }  }

  genetic.coords <- genetic.coords[-which(duplicated(sites)),]
  genetic.coords.names <- genetic.coords.names[-which(duplicated(sites))]
  differentiation <- differentiation.i
  differentiation2 <- differentiation.i2
  differentiation3 <- differentiation.i3
  
}

diag(differentiation) <- NA
diag(differentiation2) <- NA
diag(differentiation3) <- NA

plot(cost.surface,col=c("white","#75C2D3"))
points(genetic.coords[c(4,27),],col="#565656",pch=19)
lines( shortestPath(Conductance.distance, as.matrix(genetic.coords[4,]) , as.matrix(genetic.coords[27,]) , output="SpatialLines") )

plot(cost.surface,col=c("white","#75C2D3"))
points(genetic.coords[c(4,27),],col="#565656",pch=19)
lines( shortestPath(Conductance.mean.active, as.matrix(genetic.coords[4,]) , as.matrix(genetic.coords[27,]) , output="SpatialLines") )

plot(cost.surface,col=c("white","#75C2D3"))
points(genetic.coords[c(4,22),],col="#565656",pch=19)
lines( shortestPath(Conductance.mean.active, as.matrix(genetic.coords[1,]) , as.matrix(genetic.coords[13,]) , output="SpatialLines") )
lines( shortestPath(Conductance.mean.active, matrix(c(-5000000,-5e+06),ncol=2) , matrix(c(100000,-5e+06),ncol=2) , output="SpatialLines") )
lines( shortestPath(Conductance.mean.passive, matrix(c(-5000000,-5e+06),ncol=2) , matrix(c(100000,-5e+06),ncol=2) , output="SpatialLines") )

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

results.diff2 <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric( differentiation2[comb.all[i,1],comb.all[i,2]] ) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.diff3 <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric( differentiation3[comb.all[i,1],comb.all[i,2]] ) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.geoDistance <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {

  res <- as.numeric(costDistance(Conductance.distance, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))
  if(res == 0) { res <- spDistsN1( as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) , longlat=FALSE) }
  return( res )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.currentsDistance.active <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.mean.active, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.currentsDistance.passive <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(Conductance.mean.passive, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.thermalDistance.active <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(ConductanceThermalA, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## ---------------

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

results.thermalDistance.passive <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {
  
  return( as.numeric(costDistance(ConductanceThermalP, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) ))  )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## -------------------------------------------

rocky <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Genetic differentiation of Laminaria digitata across Brittany/Data/Rocky.tif"
rocky <- raster(rocky)

bathymetryHD <- raster("/Volumes/Bathyscaphe/Bio-ORACLE/Spatial Data/Netcdf/GEBCO_2019.nc")
bathymetryHD <- crop(bathymetryHD,rocky)
bathymetryHD <- resample(bathymetryHD,rocky)

bathymetryHDReclass <- bathymetryHD
bathymetryHDReclass[bathymetryHDReclass < 2 ] <- NA
bathymetryHDReclass[bathymetryHDReclass >= 2 ] <- 1

rocky <- mask(rocky,bathymetryHDReclass)

cl.2 <- makeCluster(number.cores)
registerDoParallel(cl.2)

habitat.continuity <- foreach(i=1:nrow(comb.all), .combine=c , .verbose=FALSE, .packages=c("gdistance","raster","data.table","reshape2")) %dopar% {

  line.r <- shortestPath(Conductance.distance, as.matrix(genetic.coords[comb.all[i,1],]) , as.matrix(genetic.coords[comb.all[i,2],]) , output="SpatialLines")
  result.t <- rocky[unlist(cellFromLine(rocky, line.r))]
  return( sum(result.t == 1,na.rm=T) /length(result.t) )
  
}

stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)

## -------------------------------------------

results <- data.frame(comb.all,
                      geoDistance=results.geoDistance,
                      thermalDistance.active=results.thermalDistance.active,
                      thermalDistance.passive=results.thermalDistance.passive,
                      currentsDistance.active=results.currentsDistance.active,
                      currentsDistance.passive=results.currentsDistance.passive,
                      habitat.continuity.2 = habitat.continuity.2,
                      habitat.continuity.5 = habitat.continuity.5,
                      habitat.continuity.10 = habitat.continuity.10,
                      diff=results.diff,
                      diff2=results.diff2,
                      diff3=results.diff3)

results[results == Inf] <- NA
head(results)

results[,1] <- genetic.coords.names[results[,1]]
results[,2] <- genetic.coords.names[results[,2]]

write.csv(results,file="../Final_Matrix.csv",quote = FALSE,row.names = FALSE)

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
                                         FST = as.numeric(as.character(results$diff[t.1])),
                                         NEI = as.numeric(as.character(results$diff2[t.1])),
                                         MIGRA = as.numeric(as.character(results$diff3[t.1])),
                                         
                                         Distance = results$geoDistance[t.1],

                                         habitatContinuity.2m = results$habitat.continuity.2[t.1],
                                         habitatContinuity.5m = results$habitat.continuity.5[t.1],
                                         habitatContinuity.10m = results$habitat.continuity.10[t.1],
                                         
                                         Connectivity.min.passive = min(c(results$currentsDistance.passive[t.1],results$currentsDistance.passive[t.2]),na.rm=T) ,
                                         Connectivity.mean.passive = mean(c(results$currentsDistance.passive[t.1],results$currentsDistance.passive[t.2]),na.rm=T) ,
                                         Connectivity.max.passive = max(c(results$currentsDistance.passive[t.1],results$currentsDistance.passive[t.2]),na.rm=T) ,
                                         
                                         thermalDistance.min.passive = min(c(results$thermalDistance.passive[t.1],results$thermalDistance.passive[t.2]),na.rm=T) ,
                                         thermalDistance.mean.passive = mean(c(results$thermalDistance.passive[t.1],results$thermalDistance.passive[t.2]),na.rm=T) ,
                                         thermalDistance.max.passive = max(c(results$thermalDistance.passive[t.1],results$thermalDistance.passive[t.2]),na.rm=T) ,
                                         
                                         Connectivity.min.active = min(c(results$currentsDistance.active[t.1],results$currentsDistance.active[t.2]),na.rm=T) ,
                                         Connectivity.mean.active = mean(c(results$currentsDistance.active[t.1],results$currentsDistance.active[t.2]),na.rm=T) ,
                                         Connectivity.max.active = max(c(results$currentsDistance.active[t.1],results$currentsDistance.active[t.2]),na.rm=T) ,
                                         
                                         thermalDistance.min.active = min(c(results$thermalDistance.active[t.1],results$thermalDistance.active[t.2]),na.rm=T) ,
                                         thermalDistance.mean.active = mean(c(results$thermalDistance.active[t.1],results$thermalDistance.active[t.2]),na.rm=T) ,
                                         thermalDistance.max.active = max(c(results$thermalDistance.active[t.1],results$thermalDistance.active[t.2]),na.rm=T) ,
                                         
                                         stringsAsFactors = FALSE ) )
  
}

## ---------------

# results.final[,3] <- results.final [,3] / (1 - results.final [,3])

results.final[results.final[,3] < 0,3] <- 0
results.final[results.final[,4] < 0,4] <- 0

results.final <- results.final[,c(1:9,16,20)]
colnames(results.final) <- c("pairFrom","pairTo","Fst","Nei","Migrate","marineDistance","habitatContinuity2","habitatContinuity5","habitatContinuity10","currentsDistance","thermalDistance")

results.final <- results.final[complete.cases(results.final),]

#results.final[,5] <- sqrt(results.final [,5])
#results.final[,6] <- sqrt(results.final[,6] )
#results.final[,7] <- sqrt(results.final[,7] )

fit.ibd <- lm(FST ~ Distance, data=results.final , na.action = na.omit)
fit.min.active <- lm(FST ~ Connectivity.min.active, data=results.final , na.action = na.omit)
fit.mean.active <- lm(FST ~ Connectivity.mean.active, data=results.final , na.action = na.omit)
fit.max.active <- lm(FST ~ Connectivity.max.active, data=results.final , na.action = na.omit)
fit.thermal.min <- lm(FST ~ thermalDistance.min.active, data=results.final , na.action = na.omit)
fit.thermal.mean <- lm(FST ~ thermalDistance.mean.active, data=results.final , na.action = na.omit)
fit.thermal.max <- lm(FST ~ thermalDistance.max.active, data=results.final , na.action = na.omit)

data.frame(aic.ibd= AIC(fit.ibd),
           r2.fit.ibd = summary(fit.ibd)$adj.r.squared,
           fit.min.active= AIC(fit.min.active),
           r2.fit.min.active = summary(fit.min.active)$adj.r.squared,
           fit.mean.active= AIC(fit.mean.active),
           r2.fit.mean.active = summary(fit.mean.active)$adj.r.squared,
           fit.max.active= AIC(fit.max.active),
           r2.fit.max.active = summary(fit.max.active)$adj.r.squared,
           cor.thermal.min= AIC(fit.thermal.min),
           r2.fit.thermal.min = summary(fit.thermal.min)$adj.r.squared,
           cor.thermal.mean= AIC(fit.thermal.mean),
           r2.fit.thermal.mean = summary(fit.thermal.mean)$adj.r.squared,
           cor.thermal.max= AIC(fit.thermal.max),
           r2.fit.thermal.max = summary(fit.thermal.max)$adj.r.squared
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

# write.csv(results.final,file="../Final_Pairs_Matrix3.csv",quote = FALSE,row.names = FALSE)

# save.image(file='myEnvironment.RData')
# load('myEnvironment.RData')

fit.mix <- lm(FST ~ habitatContinuity.5m+thermalDistance.mean.active+Connectivity.mean.active+Distance, data=results.final)
st <- step(fit.mix, direction = c("both"))
summary(st)

fit.final <- lm(FST ~ thermalDistance.mean.active+habitatContinuity.5m+Connectivity.mean.active, data=results.final)
AIC(fit.final)
summary(fit.final)
cor(results.final$FST ,predict(fit.final))
plot(results.final$FST ,predict(fit.final))

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