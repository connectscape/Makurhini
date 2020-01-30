#' Climate stability index

#' @param stable RasterLayer. Stable Holdridge Life Zones.
#' @param change RasterLayer. Magnitude of change in Holdridge Life Zones.
#' @param plot Logical. If TRUE, Climate stability index is plotted.
#' @param reclass Vector. Classification system where the maximum value of magnitude of change in Holdridge Life Zones is equal to 0 (Optional). Use it in case the maximum value is greater than 8, e.g., if maximum value = 9, then: reclass=c(-Inf,0,NA, 0,1,9, 1,2,8, 2,3,7, 3,4,6, 4,5,5, 5,6,4, 6,7,3, 7,8,2, 8,Inf,1).
#' @param write Character. Write the output raster e.g., "C:/folder/raster.tif"
#' @return Raster with values from 0 to 1, where 1 represents a zone of stable climate and decreases to 0 as the climate change is greater.
#' @export
CSIndex<-function(stable, change, plot=FALSE, reclass=NULL, write=NULL){
  Max_ZC <- raster::maxValue(change)
  ZCE_REC <- raster::reclassify(stable,c(-Inf,1,NA, 1,Inf,Max_ZC+1))
  if (is.null(reclass)){
    ZC_REC <- raster::reclassify(change, reclass)
    } else {
      rec <- (paste(cbind(seq(0,(Max_ZC-1), by=1)),",",cbind(c(seq(1,(Max_ZC-1), by=1), Inf)),",",cbind(seq(Max_ZC, 1)),","))
      rec <- unlist(strsplit(e, ",", fixed = TRUE))
      rec <- c(-Inf,0,NA,e)
      rec <- as.vector(as.numeric(e))
      ZC_REC <- raster::reclassify(change,rec)
      }
  CSI <- raster::mosaic(ZC_REC, ZCE_REC, fun=min)
  CSI01 <- sdmvspecies::rescale(CSI)
  if (is.null(write)){
    if(isTRUE(plot)){
      color <- RColorBrewer::colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))(100)
      raster::plot(CSI01,col = color, main = "Climate stability index", interpolate = TRUE)
      }
    } else {
      raster::writeRaster(CSI01, write, overwrite = TRUE)
      if(isTRUE(plot)){
        color <- RColorBrewer::colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))(100)
        raster::plot(CSI01, col = color, main = "Climate stability index", interpolate = TRUE)
      }
      }
  return(CSI01)
  }

