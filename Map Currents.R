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

number.cores <- 2

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm"
# directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Past LGM/Surface/LongTerm"

## --------------------------------------

xmin <- -175 ; xmax <- -104.75 ; ymin <- 20  ; ymax <- 62.5
xmin <- -180 ; xmax <- 180 ; ymin <- -70  ; ymax <- 0
xmin <- -180 ; xmax <- 180 ; ymin <- -75  ; ymax <- 90

resolution <- 0.05

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

u <- raster(Meridional.Current[4])
v <- raster(Zonal.Current[4])

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

w2hd.1 <- raster::resample(u,region.as.raster, method="bilinear")
w2hd.2 <- raster::resample(v,region.as.raster, method="bilinear")
w2hd <- brick(w2hd.1, w2hd.2)

## ----------

w2hd.1.t <- crop(w2hd.1,extent(c(-180,-30,-90,90)))
extent(w2hd.1.t) <- extent(w2hd.1.t) + c(210,210,0,0)
w2hd.1.t2 <- crop(w2hd.1,extent(c(-30,180,-90,90)))
extent(w2hd.1.t2) <- extent(w2hd.1.t2) - c(150,150,0,0)
w2hd.1 <- mosaic(w2hd.1.t,w2hd.1.t2,fun=mean,na.rm=T)

w2hd.2.t <- crop(w2hd.2,extent(c(-180,-30,-90,90)))
extent(w2hd.2.t) <- extent(w2hd.2.t) + c(210,210,0,0)
w2hd.2.t2 <- crop(w2hd.2,extent(c(-30,180,-90,90)))
extent(w2hd.2.t2) <- extent(w2hd.2.t2) - c(150,150,0,0)
w2hd.2 <- mosaic(w2hd.2.t,w2hd.2.t2,fun=mean,na.rm=T)

w2hd <- brick(w2hd.1, w2hd.2)

## ----------

# plot(w2hd)
# vectorplot(w2hd * 10, isField = "dXY", region = FALSE, margin = FALSE, narrows = 10000)

slope <- sqrt((w2hd[[1]]*10)^2 + (w2hd[[2]]*10)^2)
aspect <- atan2(w2hd[[1]]*10, w2hd[[2]]*10)
vectorplot(w2hd * 10, isField = "dXY", region = slope/10, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 10)), 
           narrows = 10000, at = seq(min(getValues(slope),na.rm=T)/10,max(getValues(slope),na.rm=T)/10,length.out=10), col.arrows='white')

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

polar = CRS("+init=epsg:3031") # North Pole polar = CRS("+init=epsg:3995")

w2hd <- raster::projectRaster(w2hd, crs = polar, method = "ngb")
slope <- raster::projectRaster(slope, crs = polar, method = "ngb")

vectorplot(w2hd, isField = "dXY", region = slope, xlab=NULL, ylab=NULL, margin = FALSE, 
           lwd.arrows=1.2, narrows = 1500 , col.arrows='white' ,  
           unit='radians', length=unit(5e-2, 'npc'),par.settings = rasterTheme(region = matlab.like(n = 50)),  at = seq(0,max(getValues(slope),na.rm=T),length.out=50)) 

## ------------------