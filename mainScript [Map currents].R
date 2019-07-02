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

number.cores <- 16

directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm"
directory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Past LGM/Surface/LongTerm"

## --------------------------------------

xmin <- 0 ; xmax <- 25 ; ymin <- -40 ; ymax <- -10

resolution <- 0.025

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

slope <- log(sqrt(w2hd[[1]]^2 + w2hd[[2]]^2)) + ( min(getValues(log(sqrt(w2hd[[1]]^2 + w2hd[[2]]^2))),na.rm=T) ) * (-1) # + 9.02 #
aspect <- atan2(w2hd[[1]], w2hd[[2]])
vectorplot(w2hd * 10, isField = "dXY", region = slope, xlab=NULL, ylab=NULL, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 100)),  narrows = 10000, at = seq(0,max(getValues(slope),na.rm=T),length.out=100)) 

## -----------------

thermal <- raster("/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Present/Surface_HR/LongTerm/Ocean.temperature.Surface.Var.Lt.Max.tif")
thermal <- crop(thermal,extent(region.as.raster))
thermal <- resample(thermal,region.as.raster,method="bilinear")

speed <- ( 1 / ( thermal / max(getValues(thermal),na.rm=TRUE)) - 1)
direction <- terrain(thermal, opt=c('aspect'), unit='degrees')

vectorplot(thermal, region = thermal, xlab=NULL, ylab=NULL, margin = FALSE, par.settings = rasterTheme(region = matlab.like(n = 100)),  narrows = 10000, at = seq(0,max(getValues(thermal),na.rm=T),length.out=100)) 
   
## -----------------

# UNDER CONSTRUCTION

direction <- (aspect / pi) * 180 # direction <- 180 + 180 * atan2(v,u) / pi
speed <- slope 

# polar = CRS("+init=epsg:3031") # North Pole polar = CRS("+init=epsg:3995")
# direction <- raster::projectRaster(direction, crs = polar, method = "ngb")
# speed <- raster::projectRaster(speed, crs = polar, method = "ngb")

direction.i <- aggregate(direction,10)
speed.i <- aggregate(speed,10)

direction.i <- direction
speed.i <- speed.i

library(geosphere)

df <- data.frame(x = as.data.frame(direction.i,xy=T)[,1], 
                 y = as.data.frame(direction.i,xy=T)[,2],
                 bearing = as.data.frame(direction.i,xy=T)[,3],
                 speed = log(as.data.frame(speed.i,xy=T)[,3]) + 9 * 10 )

df <- cbind(df,destPoint(p = df[,1:2], b = df[,3], d = df[,4]))
df <- df[complete.cases(df),]
df <- df[round(seq(1,nrow(df),length.out=10000)),]

ggplot() + geom_segment(data=df, mapping=aes(x=x, y=y, xend=lon, yend=lat), arrow=arrow(length = unit(0.2,"cm")))

ggplot(df, 
       aes(x = x , 
           y = y, 
           angle = bearing, 
           radius = speed)) +
  #geom_tile(data = gplot_data(speed.i), aes(x = x, y = y, fill = value)) +
  geom_spoke(data = df,    # this is the only difference in the plotting code
             arrow = arrow(length = unit(.05, 'inches'))) + 
  coord_equal(expand = 0) + 
  theme(legend.position = 'bottom', 
        legend.direction = 'horizontal')
  
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

polar = CRS("+init=epsg:3031") # North Pole polar = CRS("+init=epsg:3995")

# slope <- raster::projectRaster(direction, crs = polar, method = "ngb")
# aspect <- raster::projectRaster(speed, crs = polar, method = "ngb")

## ------------------
