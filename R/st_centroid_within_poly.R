#'st_centroid_within_poly
#'
#' Polygon centroid
#' @param poly object of class sf, sfc, sfg or SpatialPolygons
#' @references
#' \url{https://stackoverflow.com/questions/52522872/r-sf-package-centroid-within-polygon}
#' \url{https://stackoverflow.com/users/3609772/mitch}.
#' @importFrom magrittr %>%
#' @importFrom sf st_centroid st_within st_point_on_surface st_as_sf
#' @return
#' @export
st_centroid_within_poly <- function(poly){
  options(warn = -1)
  cl <- class(poly)[1]
  if(cl != "sf"){
    poly <- st_as_sf(poly)
  }
  centroid <-st_centroid(poly, of_largest_polygon = TRUE)
  in_poly <- st_within(centroid, poly, sparse = FALSE)
  in_poly <- diag(in_poly)
  #
  if(length(unique(in_poly)) > 1){
    st_geometry(centroid[which(in_poly == FALSE),]) <- st_geometry(st_point_on_surface(poly[which(in_poly == FALSE),]))
  }
  #
  return(centroid)
  }
