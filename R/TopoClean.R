#' Topological error correction
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param xsimplify logical or numeric.
#' @importFrom sf st_as_sf st_is_valid st_is_empty st_zm
#' @importFrom terra buffer vect
#' @importFrom rmapshaper ms_simplify
#' @export

TopoClean <- function(x, xsimplify = FALSE) {
  if(class(x)[1] == "SpatialPolygonsDataFrame") {
    x <- st_as_sf(x)
  }

  x2 <- st_zm(x, drop = TRUE)

  if(isTRUE(xsimplify) | is.numeric(xsimplify)){
    x2 <- ms_simplify(x2, keep = if(isTRUE(xsimplify)){0.9} else {xsimplify},
                      method = "vis", keep_shapes = TRUE)
  }

  x2 <- terra::buffer(vect(x2), 0)
  x2 <- st_as_sf(x2)
  x2 <- x2[which(st_is_valid(x2)),]
  x2 <- x2[which(!st_is_empty(x2)),]

  return(x2)
}
