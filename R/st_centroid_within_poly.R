#'st_centroid_within_poly
#'
#'Polygon centroid
#' @param poly object of class sf, sfc, sfg or SpatialPolygons
#' @references
#' \url{https://stackoverflow.com/questions/52522872/r-sf-package-centroid-within-polygon}
#' \url{https://stackoverflow.com/users/3609772/mitch}.
#' @importFrom magrittr %>%
#' @importFrom sf st_centroid st_within st_point_on_surface
#' @return
#' @export

st_centroid_within_poly <- function(poly){
  centroid <- poly %>% st_centroid()
  in_poly <- st_within(centroid, poly, sparse = FALSE)[[1]]
  if (in_poly){
    return(centroid)}
  centroid_in_poly <- st_point_on_surface(poly)
  return(centroid_in_poly)
  }
