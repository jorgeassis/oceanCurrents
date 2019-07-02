
require(utils)
require(colorRamps)
require(ncdf4)
require(raster)
require(rasterVis) 
library(rWind)
library(gdistance)
library(doParallel)

relocate.coordinates.na <- function(coords,rasters,maximum.distance) {
  
  set.seed(42)
  
  to.relocate <- which(raster::extract(rasters,coords[,1:2]) == 0 , arr.ind = TRUE)
  coordinates.to.relocate <- coords[to.relocate,]
  correct.points <- as.data.frame(subset(rasters,1),xy=TRUE)
  correct.points <- correct.points[correct.points[,3] !=0,1:2]
  
  if( nrow(coordinates.to.relocate) > 0 ) { 
    
    cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
    cat( paste0("\n"))
    
    near.cells <- numeric(nrow(coordinates.to.relocate))
    
    for(p in 1:nrow(coordinates.to.relocate)) {
      
      near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,1:2]),longlat=TRUE)
      
      if( near.cell.p[which.min(near.cell.p)] <= maximum.distance ) {
        
        near.cell.p <- which.min(near.cell.p)
        
      } else {   near.cell.p <- NA }
      
      
      near.cells[p] <- near.cell.p
      
    }
    
    relocated <- which(!is.na(near.cells))
    
    if( length(relocated) > 0) {
      
      near.cells <- correct.points[near.cells[relocated],]
      coords[to.relocate[relocated],] <- near.cells
      
      colnames(coords) <- c("Lon","Lat")

    }
    
  }
  
  if( nrow(coordinates.to.relocate) == 0) { 
    
    cat( paste0("None to Relocate"))
    cat( paste0("\n"))
  }
  
  to.remove <- which(raster::extract(rasters,coords[,1:2]) == 0, arr.ind = TRUE)
  
  if( length(to.remove) > 0) { 
    
    coords <- coords[-to.remove,]
    
  }
  
  ## -----------------------
  
  return( coords )
  
}




