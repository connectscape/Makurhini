#'intersect polygons
#'
#' @param p1 object of class sf, sfc, sfg or SpatialPolygons. Polygons 1 (polygons to be parallelized)
#' @param p2 object of class sf, sfc, sfg or SpatialPolygons. Polygons 2
#' @export
#' @import parallel
PIntercsect<-function(p1,p2){
  if(class(p1)[1]=="sf"){
    p1=as(p1, 'Spatial')
  } else if  (class(p2)[1]=="sf"){
    p2=as(p2, 'Spatial')
  }
  p1$idn<-rep(seq(1:round(nrow(p1)/3,0)),times=5, len=nrow(p1))
  #
  #p1$idn<-sample.int(length(p1)/2,length(p1),replace=T)
  listIZ<-split(p1, p1@data[,ncol(p1)])
  #
  cl<-parallel::makeCluster(detectCores()-1)
  listIZ<-parallel::parLapply(cl, listIZ, raster::intersect, y=p2)
  parallel::stopCluster(cl)
  p3<- do.call(rbind,listIZ)
  p3
}


