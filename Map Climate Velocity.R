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

number.cores <- 8

xmin <- -6.5 ; xmax <- 2.6 ; ymin <- 46.32  ; ymax <- 51.9
resolution <- 0.01

## -----------------

region.as.table <- matrix( NA ,nrow= ((ymax-ymin)/resolution) ,ncol= ((xmax-xmin)/resolution) )
region.as.raster <- raster(region.as.table)
extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)
crs(region.as.raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

thermal <- raster("/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm/Ocean.temperature.Surface.Var.Lt.Max.tif")
thermal <- crop(thermal,extent(region.as.raster))
thermal <- resample(thermal,region.as.raster,method="bilinear")

speed <- thermal - min(getValues(thermal),na.rm=TRUE)
speed <- speed / max(getValues(speed),na.rm=TRUE)

vectorplot(speed, region = speed, xlab=NULL, ylab=NULL, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 100)),  narrows = 5000, at = seq(min(getValues(speed),na.rm=T),max(getValues(speed),na.rm=T),length.out=100)) 

## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------