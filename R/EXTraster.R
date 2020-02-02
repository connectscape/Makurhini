#'Function to have two raster with the same number of rows, columns and extension
#' @param R.Emax RasterLayer. Raster with the largest extension. The extent of this raster will be transformed.
#' @param R.Emin RasterLayer. Raster with the least extension.
#' @export
#' @importFrom raster crop resample
EXTraster<-function(R.Emax,R.Emin){
  a <- raster::crop(R.Emax,R.Emin)
  b <- raster::resample(a,R.Emin)
  return(b)
}
