#' Topological error correction
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @importFrom sf st_as_sf st_is_valid st_is_empty st_zm st_buffer
#' @export

TopoClean <- function(x, xsimplify = FALSE) {
  if(class(x)[1] == "SpatialPolygonsDataFrame") {
    x <- st_as_sf(x)
  }

  x2 <- st_zm(x, drop = TRUE)

  if(isTRUE(xsimplify)){
    x2 <- ms_simplify(x2, 0.9, method = "vis", keep_shapes = TRUE)
  }

  x2 <- st_buffer(x2, 0)
  x2 <- x2[which(st_is_valid(x2)),]
  x2 <- x2[which(!st_is_empty(x2)),]

  return(x2)
}
