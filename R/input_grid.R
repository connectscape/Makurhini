#' ProtConn: input files
#'
#' @param node object of class sf, sfc, sfg or SpatialPolygons
#' @param landscape object of class sf, sfc, sfg or SpatialPolygons
#' @param unit character
#' @param bdist numeric
#' @param xsimplify logical
#' @importFrom sf st_sf st_cast st_buffer st_difference st_area st_geometry st_zm
#' @importFrom magrittr %>%
#' @importFrom rmapshaper ms_dissolve ms_simplify ms_clip
#' @import methods
input_grid <- function(node, landscape = NULL, unit = "ha", bdist = NULL, xsimplify = FALSE){
  class_cache <- new.env(parent = emptyenv())
  MK_Pol <- setClass("input_grid", slots = list(nodes = "sf", region = "sf",
                                                         area_unit = "character"),
                              where = class_cache)

  if(class(landscape)[1] == "SpatialPolygonsDataFrame" | class(landscape)[1] == "sf"){
    landscape <- TopoClean(landscape, xsimplify = xsimplify)
  } else{
    stop("landscape class should be SpatialPolygonDataframe or sf")
  }

  #
  mask1 <- ms_simplify(landscape, method = "vis", keep_shapes = TRUE)%>% st_buffer(., bdist)
  node <- over_poly(node, mask1, geometry = TRUE)


  if (class(node)[1] == "SpatialPolygonsDataFrame" | class(node)[1] == "sf"){
    node <- TopoClean(node, xsimplify = xsimplify)
    node <- st_cast(node, "POLYGON")
    node$IdTemp <- 1:nrow(node)
    node <- node[,which(names(node) != "geometry")]
  } else {
    stop("node class should be SpatialPolygonDataframe or sf")
  }

  if(is.null(unit)){
    unit = "ha"
  }

  MK_result <- new("input_grid", nodes = node, region = landscape, area_unit = unit)

  return(MK_result)
}



