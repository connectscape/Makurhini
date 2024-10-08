#' over sf
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param y Object of class sf, sfc, sfg or SpatialPolygons
#' @param geometry logical
#' @importFrom sf st_as_sf st_intersects st_geometry
#' @keywords internal

over_poly <- function(x, y, geometry = FALSE) {
  if(class(x)[1] == "SpatialPolygonsDataFrame") {
    x <- st_as_sf(x)
  }
  if(class(y)[1] == "SpatialPolygonsDataFrame") {
    y <- st_as_sf(y)
  }
  x.1 <- sapply(st_intersects(x, st_geometry(y)), function(z) if (length(z)==0) NA_integer_ else z[1])
  if(isTRUE(geometry)){
    x <- x[which(x.1 >= 1),]
    return(x)
  } else {
    return(x.1)
  }
}
