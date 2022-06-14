#' ProtConn: input files
#'
#' @param node object of class sf, sfc, sfg or SpatialPolygons
#' @param landscape object of class sf, sfc, sfg or SpatialPolygons
#' @param unit character
#' @param bdist numeric
#' @param xsimplify logical or numeric
#' @importFrom sf st_buffer
#' @importFrom rmapshaper ms_simplify
#' @importFrom methods setClass new
#' @keywords internal
input_grid <- function(node, landscape = NULL, unit = "ha", bdist = NULL, xsimplify = FALSE){
  class_cache <- new.env(parent = emptyenv())
  setClass("input_grid", slots = list(nodes = "sf", region = "sf",
                                                         area_unit = "character"),
                              where = class_cache)

  if(class(landscape)[1] == "SpatialPolygonsDataFrame" | class(landscape)[1] == "sf"){
    landscape <- TopoClean(landscape, xsimplify = xsimplify)
  } else {
    stop("landscape class should be SpatialPolygonDataframe or sf")
  }

  if(st_crs(node) != st_crs(landscape)){
    stop("nodes and landscape have different crs, remember that the files must have a projected coordinate system")
  }

  #
  mask1 <- rmapshaper::ms_simplify(landscape, method = "vis", keep_shapes = TRUE)
  mask1 <- st_buffer(mask1, bdist)
  node <- over_poly(node, mask1, geometry = TRUE)

  if (class(node)[1] == "SpatialPolygonsDataFrame" | class(node)[1] == "sf"){
    if(nrow(node)>0){
      node <- TopoClean(node, xsimplify = xsimplify)
      node$IdTemp <- 1:nrow(node)
      node <- node[,which(names(node) != "geometry")]
    }
  } else {
    stop("node class should be SpatialPolygonDataframe or sf")
  }

  if(is.null(unit)){
    unit = "ha"
  }

  MK_result <- new("input_grid", nodes = node, region = landscape, area_unit = unit)

  return(MK_result)
}
