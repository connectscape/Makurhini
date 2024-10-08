#' Number of vertices
#'
#' @param poly Object of class sf, sfc, sfg or SpatialPolygons
#' @importFrom sf st_cast st_as_sf
#' @keywords internal
num_vert <- function(poly){
  if(class(poly)[1] == "SpatialPolygonsDataFrame") {
    poly <- st_as_sf(poly)
  }
  poly <- st_cast(poly$geometry, "MULTIPOINT")
  poly <- sapply(poly, length)
  poly <- sum(poly)/2
  return(poly)
}
