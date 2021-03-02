#' Euclidian distances
#'
#' @param x object of class sf, sfc, sfg or SpatialPolygons. The shapefile must be in a projected coordinate system.
#' @param id character. Column name with the core id.
#' @param type_distance character. Choose one of the distances: "centroid" (faster, default), where Euclidean distance is calculated from feature centroid; "edge", where Euclidean distance is calculated from feature edges;
#' @param distance_unit character. Set a distance unit, "Makurhini::unit_covert()" compatible unit ("m", "km", "inch", "foot", "yard", "mile"). Default equal to meters "m".
#' @param keep numeric. Argument for higher processing speed. In case you have selected the "edge" distance, use this option to simplify the geometry and reduce the
#'  number of vertices (from rmapshaper::ms_simplif). The value can range from 0 to 1 and is the proportion of points to retain (default 0.02). The higher the value,
#'   the higher the speed but the greater uncertainty.
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param edgeParallel logical. Parallelize the edge distance using furrr package and multiprocess plan, default = FALSE.
#' @param pairwise logical. If TRUE a pairwise table is returned (From, To, distance) otherwise it will be a matrix.
#' @param write_table character. "" indicates output to the console.
#' @return Pairwise Euclidean distance table
#' @details The function builds on functions out of Roger Bivand and collaborators ’rgeos’ package.
#' @references Douglas, David H. and Peucker, Thomas K. (1973) "Algorithms for the Reduction of the Number of Points Required to Represent a Digitized Line or its Caricature", The Canadian Cartographer, 10(2), pp112-122.
#' @importFrom rgeos gCentroid gDistance
#' @importFrom rmapshaper ms_simplify
#' @importFrom methods as
#' @importFrom utils combn write.table
#' @importFrom future multiprocess plan availableCores
#' @importFrom furrr future_map

euclidean_distances <- function(x, id, type_distance = "centroid", distance_unit = "m",
                                keep = NULL, threshold = NULL, edgeParallel = FALSE,
                                pairwise = TRUE, write_table = NULL){
  if(missing(id)){
    stop("missing id")
  }
  if(missing(x)){
    stop("missing x object")
  }

  if(class(x)[1] == "sf") {
    x <- as(x, 'Spatial')
  }

  if(!is.null(write_table)) {
    if(!dir.exists(dirname(write_table))) {
      stop("error, output folder does not exist")
    }
  }

  if (type_distance ==  "centroid"){
    centroid_1 <- gCentroid(x, byid = TRUE)
    distance <- gDistance(centroid_1, byid = TRUE)

  } else if (type_distance == "edge"){
    if(isFALSE(edgeParallel)){
      if(!is.null(keep)){
        x_id <- x@data[,which(colnames(x@data) == id)]
        x <- ms_simplify(input = x, keep = keep, keep_shapes = TRUE, explode = FALSE)
        x$id <- x_id
        names(x) <- id
      }
      distance <- gDistance(x, byid = TRUE)
    } else {
      i = 0
      j = 0
      ng = round(nrow(x)/as.numeric(availableCores())-1)
      x2 <- list()

      repeat {
        i <- i + round(ng)
        ii <- i-(ng-1)
        if(i > nrow(x)){
          i <- nrow(x)
        }
        r <- x[ii:i,]
        j <- j+1
        x2[[j]] <- r

        if (i == nrow(x)){
          break
        }
      }
      rm(r,i,ii,j)

      works <- as.numeric(availableCores())-1
      plan(strategy = multiprocess, gc = TRUE, workers = works)
      distance <- future_map(x2, function(d){
        if(!is.null(keep)){
          d <- ms_simplify(input = d, keep = keep, keep_shapes = TRUE, explode = FALSE)
        }
        distance2 <- gDistance(d, x, byid = TRUE)
        return(distance2)
      })
      close_multiprocess(works)
      distance <- do.call(cbind, distance)
    }
  } else if (type_distance == "hausdorff-edge"){
    distance <- gDistance(x, byid = TRUE, hausdorff = TRUE)
  } else {
    stop("Error, you have to choose a type_distance option")
  }

  if (distance_unit != "m"){
    distance <- unit_convert(data_unit = distance, unit_1 = "m", unit_2 = distance_unit)
  }

  name <- c(unique(x@data[,which(colnames(x@data) == id)]), "income")
  name <- name[1:(length(name)-1)]
  colnames(distance) <- name
  rownames(distance) <- name

  if(isTRUE(pairwise)){
    xy <- t(combn(colnames(distance), 2))
    distance_2 <- data.frame(xy, distance_2 = distance[xy])
    distance_2[,1] <- as.character(distance_2[,1])
    distance_2[,2] <- as.character(distance_2[,2])
    distance_2[,1] <- as.numeric(distance_2[,1])
    distance_2[,2] <- as.numeric(distance_2[,2])
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

  return(distance_2)
}
