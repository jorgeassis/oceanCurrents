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

thermal <- raster("/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm/Ocean.temperature.Surface.Var.Lt.Max.tif")
thermal <- crop(thermal,extent(region.as.raster))
thermal <- resample(thermal,region.as.raster,method="bilinear")

speed <- thermal - min(getValues(thermal),na.rm=TRUE)
speed <- speed / max(getValues(speed),na.rm=TRUE)

vectorplot(speed, region = speed, xlab=NULL, ylab=NULL, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 100)),  narrows = 10000, at = seq(0,max(getValues(speed),na.rm=T),length.out=100)) 

## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------