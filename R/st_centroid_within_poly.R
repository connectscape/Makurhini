#'st_centroid_within_poly
#'
#'Polygon centroid
#' @param poly object of class sf, sfc, sfg or SpatialPolygons
#' @references
#' \url{https://stackoverflow.com/questions/52522872/r-sf-package-centroid-within-polygon}
#' \url{https://stackoverflow.com/users/3609772/mitch}.
#' @importFrom magrittr %>%
#' @import sf

st_centroid_within_poly <- function(poly){
  centroid <- poly %>% st_centroid()
  in_poly <- st_within(centroid, poly, sparse = F)[[1]]
  if (in_poly) return(centroid)
  centroid_in_poly <- st_point_on_surface(poly)
  return(centroid_in_poly)
  }

#' Number of vertices
#'
#' @param poly Object of class sf, sfc, sfg or SpatialPolygons
#' @import sf

num_vert <- function(poly){
  if(class(poly)[1] == "SpatialPolygonsDataFrame") {
    poly <- st_as_sf(poly)
  }
  poly <- st_cast(poly$geometry, "MULTIPOINT")
  poly <- sapply(poly, length)
  poly <- sum(poly)/2
  return(poly)
}

#' over sf
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param y Object of class sf, sfc, sfg or SpatialPolygons
#' @import sf
#' @importFrom magrittr %>%

over_poly <- function(x, y) {
  if(class(x)[1] == "SpatialPolygonsDataFrame") {
    x <- st_as_sf(x) %>% st_cast("POLYGON")
  }

  if(class(y)[1] == "SpatialPolygonsDataFrame") {
    y <- st_as_sf(y) %>% st_cast("POLYGON")
  }
  over_result <- sapply(st_intersects(x %>% st_cast("POLYGON"), st_geometry(y)), function(z) if (length(z)==0) NA_integer_ else z[1])
  return(over_result)
}
