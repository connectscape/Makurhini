#' Generates metric class
#'
#' @param metric character
#' @param attribute character
#' @param thintersect numeric
#' @param distance_threshold numeric
#' @param probability numeric
#' @param transboundary numeric
#' @param distance list
#' @export
#' @importFrom sf st_sf st_cast st_buffer st_difference st_area st_geometry
#' @importFrom magrittr %>%
#' @importFrom rmapshaper ms_dissolve ms_simplify ms_clip
#' @import methods

metric_class <- function(metric = NULL,
                         attribute = NULL,
                         thintersect = NULL,
                         distance_threshold = NULL,
                         probability = NULL,
                         transboundary = NULL,
                         distance = NULL){

  class_cache <- base::new.env(parent = emptyenv())

  if(metric == "ProtConn"){
    if(is.null(attribute)){
      attribute = "Intersected area"
    }
    if(is.null(thintersect)){
      thintersect = 0
    }
    if(is.character(thintersect)){
      stop("thintersect must be NULL or numeric")
    }
    metr <- setClass("MK_Metric", slots = list(metric = "character",
                                               attribute = "character",
                                               thintersect = "numeric",
                                               distance_threshold = "numeric",
                                               probability = "numeric",
                                               transboundary = "numeric",
                                               distance = "list"),
                     where = class_cache)



    if(is.null(transboundary)){
      tr <- 100
    } else {
      tr <- transboundary
      tr[which(tr == 0)] <- 100
      }

    metr.2 <- new("MK_Metric", metric = metric,
                  attribute = attribute,
                  thintersect = thintersect,
                  distance_threshold = distance_threshold,
                  probability = if(is.null(probability)){2} else {probability},
                  transboundary = tr,
                  distance = distance)


  } else if(metric == "PC"){
    metr <- setClass("MK_Metric", slots = list(metric = "character",
                                               distance_threshold = "numeric",
                                               probability = "numeric",
                                               distance = "list"),
                     where = class_cache)
    metr.2 <- new("MK_Metric",metric = metric,
                  distance_threshold = distance_threshold,
                  probability = if(is.null(probability)){2} else {probability},
                  distance = distance)

  } else if(metric == "IIC"){
    metr <- setClass("MK_Metric", slots = list(metric = "character",
                                               distance_threshold = "numeric",
                                               distance = "list"),
                     where = class_cache)
    metr.2 <- new("MK_Metric",metric = metric,
                  distance_threshold = distance_threshold,
                  distance = distance)
  } else {
    stop("error select between: ProtConn, PC or IIC")
  }

  return(metr.2)

}

