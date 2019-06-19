## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------
##
## C.M.I.P. Data Processing 2.0
## Muffins 'n' Code
## https://github.com/jorgeassis
##
## Python dependencies 
## https://github.com/clstoulouse/motu-client-python#UsagePIP
##
##
## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------

source("Dependencies/Functions.R")

gc(reset=TRUE)

n.clusters <- 4

## --------------------------------------

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm"

Meridional.Current <- list.files(directory,pattern="Sea.water.X.velocity",full.names = TRUE)
Meridional.Current <- Meridional.Current[which(!grepl("Range", Meridional.Current))]
Meridional.Current <- Meridional.Current[which(!grepl("Summary", Meridional.Current))]

Zonal.Current <- list.files(directory,pattern="Sea.water.Y.velocity",full.names = TRUE)
Zonal.Current <- Zonal.Current[which(!grepl("Range", Zonal.Current))]
Zonal.Current <- Zonal.Current[which(!grepl("Summary", Zonal.Current))]

length(Meridional.Current) == length(Zonal.Current)

for(f in 1:length(Zonal.Current)) {
  
  Meridional.Current.t <- raster(Meridional.Current[f])
  Zonal.Current.t <- raster(Zonal.Current[f])
  
  currents <- sqrt( (Meridional.Current.t^2) + (Zonal.Current.t^2) )
  crs(currents) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  writeRaster(currents,filename=gsub(".X.",".",Meridional.Current[f]),format="GTiff",overwrite=T)
  
}

## --------------------------------------

# Produce Range

Current.Velocity.min <- list.files(paste0(directory),pattern="Sea.water.velocity",full.names = TRUE)
Current.Velocity.min <- Current.Velocity.min[grepl("Min.tif", Current.Velocity.min)]
Current.Velocity.min <- Current.Velocity.min[!grepl("Lt", Current.Velocity.min)]

Current.Velocity.max <- list.files(paste0(directory),pattern="Sea.water.velocity",full.names = TRUE)
Current.Velocity.max <- Current.Velocity.max[grepl("Max.tif", Current.Velocity.max)]
Current.Velocity.max <- Current.Velocity.max[!grepl("Lt", Current.Velocity.max)]

for(i in 1:length(Current.Velocity.min)) {
  
  min <- raster(Current.Velocity.min[i]) ; names(min)
  max <- raster(Current.Velocity.max[i]) ; names(max)
  
  name <- unlist(gregexpr("[.]",names(min)))
  
  Current.Velocity.abs.max <- calc(stack(max,min), function(x) { max(x[1],x[2]) } )
  Current.Velocity.abs.min <- calc(stack(max,min), function(x) { min(x[1],x[2]) } )
  
  Current.Velocity.range <- Current.Velocity.abs.max - Current.Velocity.abs.min
  crs(Current.Velocity.range) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  writeRaster(Current.Velocity.range,filename= gsub("Var.Min.tif","Var.Range.tif",Current.Velocity.min[i]) ,format="GTiff",overwrite=T)
  
}

## --------------------------------------

# Plot

require(utils)
require(colorRamps)
require(ncdf4)
require(raster)
require(rasterVis) 

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm"

Meridional.Current <- list.files(directory,pattern="Sea.water.X.velocity",full.names = TRUE)
Meridional.Current <- Meridional.Current[which(!grepl("Range", Meridional.Current))]
Meridional.Current <- Meridional.Current[which(!grepl("Summary", Meridional.Current))]

Zonal.Current <- list.files(directory,pattern="Sea.water.Y.velocity",full.names = TRUE)
Zonal.Current <- Zonal.Current[which(!grepl("Range", Zonal.Current))]
Zonal.Current <- Zonal.Current[which(!grepl("Summary", Zonal.Current))]

u9 <- raster(Meridional.Current[4])
v9 <- raster(Zonal.Current[4])

w <- brick(u9, v9)
projection(w) <- CRS("+init=epsg:4326")

xmin <- -180 ; xmax <- 180 ; ymax <- 0 ; ymin <- -80
resolution <- 0.1

region.as.table <- matrix( NA ,nrow= ((ymax-ymin)/resolution) ,ncol= ((xmax-xmin)/resolution) )
region.as.raster <- raster(region.as.table)
extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)
crs(region.as.raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

w2hd.1 <- raster::resample(w[[1]],region.as.raster, method="bilinear")
w2hd.2 <- raster::resample(w[[2]],region.as.raster, method="bilinear")

w2hd <- brick(w2hd.1, w2hd.2)

plot(w2hd)
vectorplot(w2hd * 10, isField = "dXY", region = FALSE, margin = FALSE, narrows = 10000)

slope <- log(sqrt(w2hd[[1]]^2 + w2hd[[2]]^2)) + 9.02
aspect <- atan2(w2hd[[1]], w2hd[[2]])
vectorplot(w2hd * 10, isField = "dXY", region = slope, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 20)),  narrows = 10000, at = seq(0,9.33,length.out=20))

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

v9 <- crop(v9,extent(-20,20,20,50))
u9 <- crop(u9,extent(-20,20,20,50))

direction <- 180 + 180 * atan2(v9,u9) / pi
direction[is.na(direction)] <- 0

speed <- sqrt(v9^2 + u9^2)
speed[is.na(speed)] <- 0

currents <- stack(direction,speed)
names(currents) <- c("wind.direction","wind.speed")
Conductance <- flow.dispersion(currents, type ="active", fun=cost.FMGS) # cost.FMGS

library(rWind)
data(wind.data)
wind <- wind2raster(wind.data)
Conductance <- flow.dispersion(wind, type="passive")





