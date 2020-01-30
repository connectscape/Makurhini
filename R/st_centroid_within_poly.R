#'st_centroid_within_poly
#'
#'Polygon centroid
#' @param poly object of class sf, sfc, sfg or SpatialPolygons
#' @references \url{https://stackoverflow.com/questions/52522872/r-sf-package-centroid-within-polygon
#' https://stackoverflow.com/users/3609772/mitch}.
#' @export
#' @importFrom magrittr %>%

st_centroid_within_poly <- function(poly){
  centroid <- poly %>% sf::st_centroid()
  in_poly <- sf::st_within(centroid, poly, sparse = F)[[1]]
  if (in_poly) return(centroid)
  centroid_in_poly <- sf::st_point_on_surface(poly)
  return(centroid_in_poly)
  }


