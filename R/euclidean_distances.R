#' Euclidian distances
#'
#' @param x object of class sf, sfc, sfg or SpatialPolygons. It must be in a projected coordinate system.
#' @param id character. Column name with the patches id.
#' @param centroid logical. If centroid is TRUE then the Euclidean distance is calculated from patch centroids, if FALSE the Euclidean distance is calculated from patch edges.
#' @param distance_unit character. Set a distance unit, "Makurhini::unit_covert()" compatible unit ("m", "km", "inch", "foot", "yard", "mile"). Default equal to meters "m".
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param keep numeric. Argument for higher processing speed. Only is available for "edge" distance or centroid equal FALSE, use this option to simplify the geometry and reduce the
#'  number of vertices. The value can range from 0 to 1 and is the proportion of points to retain (e.g., 0.02). The higher the value,
#'   the higher the speed but the greater uncertainty.
#' @param distParallel logical or numeric. Only is available for "edge" and "centroid" distance, use this option to parallelize using furrr package and multiprocess plan, default = FALSE.
#' @param ActiveParallel logical. It should be TRUE if there is already an open parallelization plan.
#' @param pairwise logical. If TRUE a pairwise table is returned (From, To, distance) otherwise it will be a matrix.
#' @param write_table character. "" indicates output to the console.
#' @return Pairwise Euclidean distance table
#' @references Douglas, David H. and Peucker, Thomas K. (1973) "Algorithms for the Reduction of the Number of Points Required to Represent a Digitized Line or its Caricature", The Canadian Cartographer, 10(2), pp112-122.
#' @importFrom sf st_as_sf st_distance st_geometry_type
#' @importFrom rmapshaper ms_simplify
#' @importFrom methods as
#' @importFrom utils combn write.table
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map
#' @export
euclidean_distances <- function(x, id, centroid = TRUE, distance_unit = "m",
                                keep = NULL, threshold = NULL, distParallel = FALSE,
                                ActiveParallel = FALSE,
                                pairwise = TRUE, write_table = NULL){
  if(missing(id)){
    stop("missing id")
  }
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(missing(x)){
    stop("missing x object")
  }

  if(class(x)[1] != "sf") {
    x <- st_as_sf(x)
  }

  if(!is.null(write_table)) {
    if(!dir.exists(dirname(write_table))) {
      stop("error, output folder does not exist")
    }
  }

  if (isTRUE(centroid)){
    if(isFALSE(distParallel) | is.null(distParallel) & isFALSE(ActiveParallel)){
      if(as.character(unique(st_geometry_type(x))) != "POINT"){
        x <- st_centroid_within_poly(x)
      }
      distance <- st_distance(x, by_element = FALSE); attr(distance, "units") <- NULL; class(distance) <- setdiff(class(distance),"units")
    } else {
      if(isFALSE(ActiveParallel)){
        if(is.logical(distParallel)){
          works <- as.numeric(availableCores())-1
        } else {
          if(is.numeric(distParallel)){
            works <- distParallel
          }
        }

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        plan(strategy = strat, gc = TRUE, workers = works)
      }

      if(as.character(unique(st_geometry_type(x))) != "POINT"){
        x <- future_map_dfr(split(x, ceiling(1:nrow(x)/10)), function(y){
          y <- st_centroid_within_poly(y); return(y)})
      }

      distance <- future_map(split(x, ceiling(1:nrow(x)/1000)), function(y){
        y <- st_distance(y, x, by_element = FALSE); attr(y, "units") <- NULL; class(y) <- setdiff(class(y),"units")
        return(y)}); distance <- do.call(rbind, distance)

      if(isFALSE(ActiveParallel)){
        close_multiprocess()
      }

    }
  } else {
    if(isFALSE(distParallel) | is.null(distParallel) & isFALSE(ActiveParallel)){
      if(!is.null(keep) & any(as.character(unique(st_geometry_type(x))) != "POINT")){
        x_id <- x[[id]]; x.1 <- tryCatch(ms_simplify(input = x, keep = keep,
                                                     keep_shapes = TRUE,
                                                     explode = FALSE),
                                         error = function(err)err)

        if(!inherits(x.1, "error")){
          if(any(st_is_empty(x.1))){
            x.1 <- rbind(x.1[which(!st_is_empty(x.1)),],
                         x[which(st_is_empty(x.1)),])
            x.1 <- x.1[sapply(x_id, function(y){which(x.1[[id]] == y)}),]
          }

          if(nrow(x.1) == nrow(x)){
            x <- x.1; x[which(names(x) == id)] <- x_id
          } else {
            x <- rbind(x.1, x[which(x.1[[id]] %!in% x[[id]])])
            x <- x[sapply(x_id, function(y){which(x[[id]] == y)}),]
            x[which(names(x) == id)] <- x_id
            message("It was not possible to simplify the shapes in the estimation of the distances between edges")
          }
        } else {
          stop("error in distParallel and keep parameter")
        }
      }
      distance <- st_distance(x, by_element = FALSE); attr(distance, "units") <- NULL; class(distance) <- setdiff(class(distance),"units")
    } else {
      if(isFALSE(ActiveParallel)){
        if(is.logical(distParallel)){
          works <- as.numeric(availableCores())-1
        } else {
          if(is.numeric(distParallel)){
            works <- distParallel
          }
        }

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        plan(strategy = strat, gc = TRUE, workers = works)
      }

      distance <- tryCatch(future_map(split(x, 1:nrow(x)), function(y){
        if(!is.null(keep) & as.character(unique(st_geometry_type(y))) != "POINT"){
          y.1 <- tryCatch(ms_simplify(input = y, keep = keep, keep_shapes = TRUE, explode = FALSE), error = function(err)err)
          if(!inherits(y.1, "error")) {
            if(any(st_is_empty(y.1))){
              y.1 <- rbind(y.1[which(!st_is_empty(y.1)),],
                           y[which(st_is_empty(y.1)),])
            }
            if(nrow(y.1) == nrow(y)){
              y <- y.1
            } else {
              y.1 <- rbind(y.1, y[which(y.1[[id]] %!in% y[[id]])])
              y <- y.1[sapply(y[[id]], function(i){which(y.1[[id]] == i)}),]
            }
          } else {
            stop("error in distParallel")
          }
        }
        distance <- st_distance(y, x, by_element = FALSE); attr(distance, "units") <- NULL; class(distance) <- setdiff(class(distance),"units")
        return(distance)
      }), error = function(err)err)

      if(isFALSE(ActiveParallel)){
        close_multiprocess()
      }

      if(inherits(distance, "error")) {
        stop(distance)
      } else {
        distance <- do.call(rbind, distance)
      }
    }
  }

  if (distance_unit != "m"){
    distance <- unit_convert(data_unit = distance, unit_1 = "m", unit_2 = distance_unit)
  }

  name <- unique(x[[which(names(x) == id)]]); colnames(distance) <- name; rownames(distance) <- name

  if(isTRUE(pairwise)){
    xy <- t(combn(colnames(distance), 2)); distance_2 <- data.frame(xy, distance_2 = distance[xy])
    distance_2[,1] <- as.character(distance_2[,1]); distance_2[,2] <- as.character(distance_2[,2])
    distance_2[,1] <- as.numeric(distance_2[,1]);distance_2[,2] <- as.numeric(distance_2[,2])
    names(distance_2) <- c("From", "To", "Distance")

    if (!is.null(threshold)){
      distance_2 <-  distance_2[which(distance_2$Distance <= threshold),]
    }
  } else {
    distance_2 <- distance
    if (!is.null(threshold)){
      distance_2[which(distance_2 <= threshold)] <- NA
    }
  }

  if(!is.null(write_table)){
    write.table(distance_2, write_table, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  invisible(gc())
  return(distance_2)
}
