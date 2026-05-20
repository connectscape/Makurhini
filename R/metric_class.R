#' Generates metric class
#'
#' @param metric character
#' @param attribute character
#' @param thintersect numeric
#' @param distance_threshold numeric
#' @param probability numeric
#' @param pij_min numeric
#' @param transboundary numeric
#' @param distance list
#' @importFrom methods setClass new
#' @keywords internal

metric_class <- function(metric = NULL,
                         attribute = NULL,
                         thintersect = NULL,
                         distance_threshold = NULL,
                         probability = NULL,
                         pij_min = NULL,
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

    if(is.null(pij_min)){
      pij_min = Inf
    }

    metr <- setClass("MK_Metric", slots = list(metric = "character",
                                               attribute = "character",
                                               thintersect = "numeric",
                                               distance_threshold = "numeric",
                                               probability = "numeric",
                                               pij_min = "numeric",
                                               transboundary = "numeric",
                                               distance = "list"),
                     where = class_cache)

    if(is.null(transboundary)){
      tr <- 100
    } else {
      tr <- transboundary; tr[which(tr == 0)] <- 100
      }

    metr <- new("MK_Metric", metric = metric,
                attribute = attribute,
                thintersect = thintersect,
                distance_threshold = distance_threshold,
                probability = if(is.null(probability)){0.5} else {probability},
                pij_min = pij_min,
                transboundary = tr,
                distance = distance)
    return(metr)
  } else if(metric == "PC"){
    if(is.null(pij_min)){
      pij_min = Inf
    }

    metr <- setClass("MK_Metric", slots = list(metric = "character",
                                               distance_threshold = "numeric",
                                               probability = "numeric",
                                               pij_min = "numeric",
                                               distance = "list"),
                     where = class_cache)
    metr <- new("MK_Metric", metric = metric,
                  distance_threshold = distance_threshold,
                  probability = if(is.null(probability)){0.5} else {probability},
                  pij_min = pij_min,
                  distance = distance)
    return(metr)
  } else if(metric == "IIC"){
    metr <- setClass("MK_Metric", slots = list(metric = "character",
                                               distance_threshold = "numeric",
                                               distance = "list"),
                     where = class_cache)
    metr <- new("MK_Metric", metric = metric,
                  distance_threshold = distance_threshold,
                  distance = distance)
    return(metr)
  } else {
    stop("error select between: ProtConn, PC or IIC")
  }
}

