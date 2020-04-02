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
#' @references Douglas, David H. and Peucker, Thomas K. (1973) "Algorithms for the Reduction of the Number of Points Required to Represent a Digitized Line or its Caricature", The Canadian Cartographer, 10(2), pp112-122.
#' @importFrom rgeos gCentroid gDistance
#' @importFrom rmapshaper ms_simplify
#' @importFrom methods as
#' @importFrom utils combn write.table

euclidean_distances <- function(x, id, type_distance = "centroid", distance_unit = "m",
                               keep = 0.05, threshold = NULL, edgeParallel = FALSE,
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
        ng = round(nrow(x)/8)
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

        plan(strategy = multiprocess, gc = TRUE)
        distance <- future_map(x2, function(d){
          if(!is.null(keep)){
          d <- ms_simplify(input = d, keep = keep, keep_shapes = TRUE, explode = FALSE)
          }
          distance2 <- gDistance(d, x, byid = TRUE)
          return(distance2)
        })
        distance <- do.call(cbind, distance)
        future:::ClusterRegistry("stop")
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

#' Cost distances
#'
#' @param x object of class sf, sfc, sfg or SpatialPolygons. The shapefile must be in a projected coordinate system.
#' @param id character. Column name with the core id.
#' @param type_distance character. Choose one of the distances: "least-cost" (default) that takes into account obstacles and local friction of the landscape (See, gdistance package);
#'  "commute-time" that is analogous to the resistance distance. This distance is based on the random walk theory and calculated using the electrical circuit theory (See, gdistance package).
#' @param resistance raster. Raster object with resistance values (landscape friction).
#' @param CostFun A function to compute the cost to move between cells.The default is the conductance
#'  (isotropic cost distance):  function(x) 1/mean(x).
#' @param ngh numeric. Neighbour graph (directions) for distance calculations: 4 (von Neu-mann neighbourhood), 8 (Moore neighbourhood) or 16 (king’s and knight’s moves). Default equal 16.
#' @param mask object of class sf, sfc, sfg or SpatialPolygons. For higher processing speed use this option to clip the resistance at the extent of the mask.
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param geometry_out numeric. Avoids the "Inf" error and corrects the "Inf" distance values when some spatial geometries are out
#' of the resistance extent. The raster NA values will be replaced by a CONDUCTANCE value that must be provided by the user
#' using this argument, so that it is possible to calculate cost distances for all the pairs of nodes.
#'  If NULL, then a Euclidean distance will be calculated to find these distances.
#' @param bounding_circles numeric. If a value is entered, this will create bounding circles around pairs of core areas
#'  (recommended for speed, large resistance rasters or pixel resolution < 150 m).
#'  Buffer distances are entered in map units. Also, the function is parallelized using
#'   and furrr package and multiprocess plan, default = NULL.
#' @param pairwise logical. If TRUE a pairwise table is returned (From, To, distance) otherwise it will be a matrix.
#' @param write_table character. "" indicates output to the console.
#' @return cost distance matrix
#' @references https://cran.r-project.org/web/packages/gdistance/gdistance.pdf
#' @importFrom magrittr %>%
#' @import sf
#' @importFrom purrr map_dbl
#' @importFrom rgeos gTouches gCentroid
#' @importFrom dismo voronoi
#' @importFrom raster crop mask extend buffer raster
#' @importFrom gdistance transition geoCorrection costDistance commuteDistance
#' @importFrom methods as new
#' @importFrom stats na.omit
#' @importFrom utils combn write.table
#' @importFrom future multiprocess plan
#' @importFrom furrr future_map

cost_distances <- function(x, id, type_distance = "least-cost", resistance = NULL, CostFun = NULL,
                           ngh = NULL, mask = NULL, threshold = NULL,
                           geometry_out = NULL,
                           bounding_circles = NULL,
                           pairwise = TRUE,
                           write_table = NULL){

  if(missing(id)){
    stop("missing id")
  }

  if(missing(x)){
    stop("missing x object")
  }

  if (!is.null(write_table)) {
    if (!dir.exists(dirname(write_table))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(x)[1] == "SpatialPolygonsDataFrame") {
    x <- st_as_sf(x)
  }

  if (is.null(resistance)){
    stop("Error, you need a resistance raster")
  }

  coordenates_1 <- x
  coordenates_1$lon <- map_dbl(x$geometry, ~ st_centroid_within_poly(.x)[[1]])
  coordenates_1$lat <- map_dbl(x$geometry, ~ st_centroid_within_poly(.x)[[2]])

  if (!is.null(mask)){
    if(class(mask)[1] == "sf"){
      mask <- as(mask, 'Spatial')
    }
    resistance_1 <- crop(resistance, mask)
  } else {
    resistance_1 <- resistance
  }

  if(is.null(ngh)){
    ngh = 16
  }

  if(is.null(CostFun)){
    CostFun <-  function(x) 1/mean(x)
  }

  if(is.null(bounding_circles)){
    st_geometry(coordenates_1) <- NULL
    coordenates_1 <- cbind(coordenates_1$lon, coordenates_1$lat)
    rownames(coordenates_1) <- NULL

   if(type_distance == "least-cost"){
     Iso_conductance <- transition(resistance_1, CostFun, directions = ngh) %>%
       geoCorrection(., scl = FALSE)
     distance_result <- costDistance(Iso_conductance, coordenates_1) #reciprocal of conductance, works with resistance
     } else if (type_distance == "commute-time") {
        Iso_conductance <- transition(resistance_1, CostFun, directions = ngh) %>%
          geoCorrection(., type = "r")
        distance_result <- commuteDistance(Iso_conductance, coordenates_1)
      } else {
        stop("Error, you have to choose a type_distance option")
      }

      distance_result <- as.matrix(distance_result)
      rownames(distance_result)<- x[[id]]
      colnames(distance_result)<- x[[id]]

      distance_result_matrix <- distance_result

      distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
      distance_result <- as.data.frame(as.table(distance_result))
      distance_result <- na.omit(distance_result)
      distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
      distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
      names(distance_result)<-c("From", "To", "Distance")

    #Inf
    x2 <- distance_result[which(is.infinite(distance_result[,3])),]

    if(nrow(x2) >= 1){
      Infinitos <- paste0(x2$From, "_", x2$To)
      distance_result2 <- distance_result[which(!is.infinite(distance_result[,3])),]

      error <- unique(x2$To)
      min_dist <- euclidean_distances(x, id, type_distance = "centroid", distance_unit = "m")
      #
      if (is.null(geometry_out)){
        distance2 <- min_dist
        distance2$idn <- paste0(distance2$From, "_", distance2$To)
        #Filter distances
        distance2 <- distance2[which(distance2$idn %in% Infinitos),]
        distance2$idn <- NULL
        names(distance_result2)[3] <- "Distance"
        distance_result <- rbind(distance_result2, distance2)


        for(i in 1:nrow(distance2)){
          distance_result_matrix[which(rownames(distance_result_matrix) == distance2$From[i]),
                                 which(colnames(distance_result_matrix) == distance2$To[i])] <- distance2[i,3]

          distance_result_matrix[which(rownames(distance_result_matrix) == distance2$To[i]),
                                 which(colnames(distance_result_matrix) == distance2$From[i])] <- distance2[i,3]
        }

      } else {
        Iso_conductance <- raster(Iso_conductance)
        Iso_conductance <- extend(Iso_conductance, as(x, 'Spatial'))

        Iso_conductance[is.na(Iso_conductance[])] <- geometry_out
        Iso_conductance <- new("TransitionLayer", Iso_conductance)

        #Get new distances
        if(type_distance == "least-cost"){
          distance_result<- costDistance(Iso_conductance, coordenates_1) #reciprocal of conductance, works with resistance
        } else if (type_distance == "commute-time") {
          distance_result <- commuteDistance(Iso_conductance, coordenates_1)
        }

        distance_result <- as.matrix(distance_result)
        rownames(distance_result)<- x[[id]]
        colnames(distance_result)<- x[[id]]

        distance_result_matrix <- distance_result

        distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
        distance_result <- as.data.frame(as.table(distance_result))
        distance_result <- na.omit(distance_result)
        distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
        distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
        names(distance_result)<-c("From", "To", "Distance")
      }
    }
  } else {
    ####voronoi vecinos

    voronoi_adj <- voronoi(gCentroid(as(x, 'Spatial'), byid = T)) %>%
      gTouches(., byid = T)
    colnames(voronoi_adj) <- x[[id]]
    rownames(voronoi_adj) <- x[[id]]

    plan(strategy = multiprocess, gc = TRUE)

    distance <- future_map(as.list(colnames(voronoi_adj)), function(i){
      ##
      foc_adj <- rbind(x[which(x[[id]] == i),],
                       x[as.vector(which(voronoi_adj[, which(x[[id]] == i)])),]) %>%
        as(., 'Spatial')

      dist <- euclidean_distances(foc_adj, id, threshold = bounding_circles, distance_unit = "m" )
      foc_adj <- foc_adj[which(foc_adj[[id]] %in% c(unique(dist$From), unique(dist$To))),]

      if(nrow(foc_adj) > 1){
        r2 <- gCentroid(foc_adj, byid = TRUE) %>% buffer(., width = bounding_circles)
        r2 <- crop(resistance_1, r2) %>% raster::mask(., r2)

        ##
        coordenates_2 <- coordenates_1[which(coordenates_1[[id]] %in% foc_adj[[id]]),]
        st_geometry(coordenates_2) <- NULL
        coordenates_2 <- cbind(coordenates_2$lon, coordenates_2$lat)
        rownames(coordenates_2) <- NULL
        ##

        if(type_distance == "least-cost"){
          Iso_conductance <- transition(r2, CostFun, directions = ngh)%>%
            geoCorrection(., scl = FALSE)
          distance_result <- costDistance(Iso_conductance, coordenates_2) #reciprocal of conductance, works with resistance
        } else if (type_distance == "commute-time") {
          Iso_conductance <- transition(r2, CostFun, directions = ngh)%>%
            geoCorrection(., type = "r")
          distance_result <- commuteDistance(Iso_conductance, coordenates_2)
        } else {
          stop("Error, you have to choose a type_distance option")
        }

        distance_result <- as.matrix(distance_result)
        rownames(distance_result)<- foc_adj[[id]]
        colnames(distance_result)<- foc_adj[[id]]

        distance_result_matrix <- distance_result

        distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
        distance_result <- as.data.frame(as.table(distance_result))
        distance_result <- na.omit(distance_result)
        distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
        distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
        names(distance_result)<-c("From", "To", "Distance")

        #Inf
        x2 <- distance_result[which(is.infinite(distance_result[,3])),]

        if(nrow(x2) >= 1){
          Infinitos <- paste0(x2$From, "_", x2$To)
          distance_result2 <- distance_result[which(!is.infinite(distance_result[,3])),]

          error <- unique(x2$To)
          min_dist <- euclidean_distances(foc_adj, id, type_distance = "centroid", distance_unit = "m")

          if (is.null(geometry_out)){
            distance2 <- min_dist
            distance2$idn <- paste0(distance2$From, "_", distance2$To)
            #Filter distances
            distance2 <- distance2[which(distance2$idn %in% Infinitos),]
            distance2$idn <- NULL
            names(distance_result2)[3] <- "Distance"
            distance_result <- rbind(distance_result2, distance2)

            for(i in 1:nrow(distance2)){
              distance_result_matrix[which(rownames(distance_result_matrix) == distance2$From[i]),
                                     which(colnames(distance_result_matrix) == distance2$To[i])] <- distance2[i,3]

              distance_result_matrix[which(rownames(distance_result_matrix) == distance2$To[i]),
                                     which(colnames(distance_result_matrix) == distance2$From[i])] <- distance2[i,3]
            }

          } else {
            Iso_conductance <- raster(Iso_conductance)
            Iso_conductance <- extend(Iso_conductance, as(x, 'Spatial'))

            Iso_conductance[is.na(Iso_conductance[])] <- geometry_out
            Iso_conductance <- new("TransitionLayer", Iso_conductance)

            #Get new distances
            if(type_distance == "least-cost"){
              distance_result<- costDistance(Iso_conductance, coordenates_2) #reciprocal of conductance, works with resistance
            } else if (type_distance == "commut|e-time") {
              distance_result <- commuteDistance(Iso_conductance, coordenates_2)
            }
            rm(Iso_conductance)

            distance_result <- as.matrix(distance_result)
            rownames(distance_result)<- foc_adj[[id]]
            colnames(distance_result)<- foc_adj[[id]]

            distance_result_matrix <- distance_result

            distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
            distance_result <- as.data.frame(as.table(distance_result))
            distance_result <- na.omit(distance_result)
            distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
            distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
            names(distance_result)<-c("From", "To", "Distance")
          }
        }
        if(isFALSE(pairwise)){
          distance_result <- distance_result_matrix
        }

        return(distance_result)
      } }, .progress = TRUE)
    future:::ClusterRegistry("stop")

    distance <- purrr::compact(distance)

    if(isFALSE(pairwise)){
      uno <- list()
      for(i in 1:length(distance)){
        dos <- rownames(distance[[i]])
        uno[[i]] <- dos
      }
      uno <- do.call(c, uno)
      uno <- unique(uno)

      mm <- matrix(nrow = length(uno), ncol = length(uno))
      rownames(mm) <- uno
      colnames(mm) <- uno

      for(i in uno){
        mm2 <- mm[,i] %>% as.matrix()

        for(j in rownames(mm2)){
          if(j != i){
            mm2[j,] <- min(map(distance, function(x){
              x2 <- x[which(rownames(x)==i),which(rownames(x)==j)]
              if(length(x2) == 0){
                x2 = NULL
              }
              return(x2)})%>% purrr::compact(.) %>% unlist())
          } else {
            mm2[j,] <- 0
          }
          mm[,i] <- mm2
        }
      }
      distance_result_matrix <- mm
    } else {
      distance_result <- do.call(rbind, distance)
      distance_result$dn <- paste0(distance_result$From, "_",distance_result$To)

      distance_result <- map(unique(distance_result$dn), function(i){
        d1 <- distance_result[which(distance_result$dn == i),]
        d1 <- d1[which(d1$Distance == min(d1$Distance)),1:3]
        if(nrow(d1)>1){
          d1 <- d1[1,]
        }
        return(d1)}) %>% do.call(rbind,.)
      rownames(distance_result) <- NULL
    }
    }


  if(isTRUE(pairwise)){
    if (!is.null(threshold)){
      distance_result <-  distance_result[which(distance_result[,3] <= threshold),]
    }
    rownames(distance_result) <- NULL
  } else {
    distance_result <- distance_result_matrix
    distance_result <- as.matrix(distance_result)
    if (!is.null(threshold)){
      distance_result[which(distance_result <= threshold)] <- NA
    }
  }

  if(!is.null(write_table)){
    write.table(distance_result, write_table, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  return(distance_result)
  }
