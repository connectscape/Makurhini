#' ProtConn: input files
#'
#' @param node object of class sf, sfc, sfg or SpatialPolygons
#' @param landscape object of class sf, sfc, sfg or SpatialPolygons
#' @param unit character
#' @param bdist numeric
#' @param xsimplify logical or numeric
#' @importFrom sf st_buffer st_crs
#' @importFrom rmapshaper ms_simplify
#' @importFrom methods setClass new
#' @keywords internal
input_grid <- function(node, landscape = NULL, unit = "ha", bdist = NULL, xsimplify = FALSE){
  class_cache <- new.env(parent = emptyenv())
  setClass("input_grid", slots = list(nodes = "sf", region = "sf",
                                                         area_unit = "character"),
                              where = class_cache)

  if(any(class(landscape)[1] == "SpatialPolygonsDataFrame" | class(landscape)[1] == "sf")){
    landscape <- TopoClean(x = landscape, xsimplify = xsimplify)
  } else {
    stop("landscape class should be SpatialPolygonDataframe or sf")
  }

  if(st_crs(node) != st_crs(landscape)){
    stop("nodes and landscape have different crs, remember that the files must have a projected coordinate system")
  }

  if(is.null(unit)){
    unit = "ha"
  }

  node <- rmapshaper::ms_simplify(landscape, method = "vis", keep_shapes = TRUE) |>
    st_buffer(x = _, bdist)|> over_poly(x = node, y = _, geometry = TRUE)

  if (class(node)[1] == "SpatialPolygonsDataFrame" | class(node)[1] == "sf"){
    if(nrow(node)>0){
      node <- TopoClean(node, xsimplify = xsimplify)
      node$IdTemp <- 1:nrow(node); node <- node[,which(names(node) != "geometry")]
      return(new("input_grid", nodes = node, region = landscape, area_unit = unit))
    } else {
      stop("No nodes, please check for possible topology errors")
    }
  } else {
    stop("Node class should be SpatialPolygonDataframe or sf")
  }
}
