#' Topological error correction
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param xsimplify logical or numeric.
#' @importFrom sf st_as_sf st_is_valid st_is_empty st_zm
#' @importFrom terra buffer vect
#' @importFrom rmapshaper ms_simplify
#' @keywords internal
TopoClean <- function(x, xsimplify = FALSE) {
  if(class(x)[1] == "SpatialPolygonsDataFrame" | class(x)[1] == "SpatVector") {
    x <- st_as_sf(x)
  }

  x <- st_zm(x, drop = TRUE)

  if(isTRUE(xsimplify) | is.numeric(xsimplify)){
    x <- ms_simplify(x, keep = if(isTRUE(xsimplify)){0.9} else {xsimplify},
                      method = "vis", keep_shapes = TRUE)
  }

  x <- terra::buffer(vect(x), 0) |> st_as_sf(x = _)
  x <- x[which(st_is_valid(x)),]; x <- x[which(!st_is_empty(x)),]
  return(x)
}
