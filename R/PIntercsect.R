#'intersect polygons
#'
#' @param p1 object of class sf, sfc, sfg or SpatialPolygons. Polygons 1 (polygons to be parallelized)
#' @param p2 object of class sf, sfc, sfg or SpatialPolygons. Polygons 2
#' @importFrom parallel makeCluster stopCluster parLapply detectCores
#' @importFrom methods as
#' @importFrom raster intersect
PIntercsect<-function(p1,p2){
  if(class(p1)[1]=="sf"){
    p1 <- as(p1, 'Spatial')
  } else if  (class(p2)[1]=="sf"){
    p2 <- as(p2, 'Spatial')
  }
  p1$idn <- rep(seq(1:round(nrow(p1)/3,0)),times=5, len=nrow(p1))
  #
   listIZ <- split(p1, p1@data[,ncol(p1)])
  #
  cl <- makeCluster(detectCores()-1)
  listIZ <- parLapply(cl, listIZ, intersect, y=p2)
  stopCluster(cl)
  p3 <- do.call(rbind,listIZ)
  return(p3)
}


