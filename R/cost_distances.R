#' Estimate isotropic cost distances
#'
#' @param x object of class sf, sfc, sfg or SpatialPolygons. The shapefile must be in a projected coordinate system (number of patches > 2).
#' @param id character. Column name with the patches id.
#' @param LCD logical. If TRUE “least-cost distance” is estimated(default) else "commute-time" , analogous to the resistance distance of circuitscape, will be estimated (See, gdistance package).
#' @param resistance raster. Raster object with resistance values (landscape friction). Consider using a range of resistance values between 1 and 10, where 10 is the maximum resistance to movement, NA should be a complete barrier.
#' @param CostFun A function to compute the cost to move between cells. The default is the conductance
#'  (isotropic cost distance): function(x) 1/x[2].
#' @param ngh numeric. Neighbor graph (directions) for distance calculations: 4 (von Neu-mann neighbourhood), 8 (Moore neighbourhood) or 16 (king’s and knight’s moves). Default equal 8.
#' @param mask object of class sf, sfc, sfg or SpatialPolygons. For higher processing speed use this option to clip the resistance at the extent of the mask.
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param geometry_out numeric. Avoids the "Inf" error and corrects the "Inf" distance values when some spatial geometries are out
#' of the resistance extent. The raster NA values will be replaced by a Resistance value that must be provided by the user
#' using this argument, so that it is possible to calculate cost distances for all the pairs of nodes.
#'  If NULL, then a Euclidean distance will be calculated to find these distances.
#' @param bounding_circles numeric. If a value is entered, this will create bounding circles around pairs of core areas
#'  (recommended for speed, large resistance rasters or pixel resolution < 150 m).
#'  Buffer distances are entered in map units. Also, the function is parallelized using
#'   and furrr package and multiprocess plan, default = NULL.
#' @param least_cost.java logical. If TRUE then the programming language and computing platform 'java' will be used
#' to estimate the least-cost distance USING only the FORMULA: function(x) 1/mean(x). It is necessary to have java installed. This option use
#' the package 'graph4lg' to reduce computation times. WARNING: This function only works with
#' integer type resistance raster (it does not accept decimals or floating point resistance).
#' @param cores.java numeric. Computer cores used to run the .jar file (see, graph4lg), default = 1.
#' @param ram.java numeric. RAM gigabytes to run the .jar file (see, graph4lg), default = NULL.
#' @param pairwise logical. If TRUE a pairwise table is returned (From, To, distance) otherwise it will be a matrix.
#' @param write_table character. "" indicates output to the console.
#' @return cost distance matrix
#' @details The function builds on functions out of Jacob van Etten’s ’gdistance’ package.
#' @references Jacob van Etten. 2017. R Package gdistance: Distances and Routes on Geographical Grids.
#' 	Journal of Statistical Software. \url{10.18637/jss.v076.i13}
#'
#' 	Paul Savary. 2020. R Package graph4lg: Build Graphs for Landscape Genetics Analysis.
#'  \url{https://cran.r-project.org/web/packages/graph4lg/index.html}
#' @importFrom sf st_as_sf st_geometry st_geometry<- st_coordinates st_centroid st_combine st_voronoi st_touches st_collection_extract
#' @importFrom rgeos gCentroid
#' @importFrom raster crop mask extend buffer raster res
#' @importFrom gdistance transition geoCorrection costDistance commuteDistance
#' @importFrom methods as new
#' @importFrom stats na.omit
#' @importFrom utils combn write.table
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map
#' @importFrom graph4lg mat_cost_dist
#' @importFrom magrittr %>%
#' @export

cost_distances <- function(x, id,
                           LCD = "least-cost",
                           resistance = NULL,
                           CostFun = NULL,
                           ngh = NULL, mask = NULL, threshold = NULL,
                           geometry_out = NULL,
                           bounding_circles = NULL,
                           least_cost.java = FALSE,
                           cores.java = 1, ram.java = NULL,
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

  if(class(x)[1] != "sf") {
    x <- st_as_sf(x)
  }

  if (is.null(resistance)){
    stop("Error, you need a resistance raster")
  }
  . = NULL
  cord <- st_centroid_within_poly(x) %>% st_coordinates(.)
  coordenates_1 <- x; coordenates_1$lon <- cord[,1]; coordenates_1$lat <- cord[,2]

  if (!is.null(mask)){
    if(class(mask)[1] == "sf"){
      mask <- as(mask, 'Spatial')
    }
    resistance_1 <- crop(resistance, mask)
  } else {
    resistance_1 <- resistance
  }

  if(is.null(ngh)){
    ngh = 8
  }

  if(is.null(CostFun)){
    CostFun <-  function(x) 1/x[2]
  }

  resistance_1@crs@projargs <- "+proj=merc +units=m"

  if(isTRUE(least_cost.java)){
    pts_1 <- data.frame(ID = x[[id]],
                        x = coordenates_1$lon,
                        y = coordenates_1$lat)

    val <- raster::values(resistance_1)
    cost <- data.frame(code = unique(val),
                       cost = unique(val))

    if(length(which(is.na(cost$cost))) >0){
      cost <- cost[-which(is.na(cost$cost)),]
    }

    distance_result <- mat_cost_dist(raster = resistance_1,
                           pts = pts_1,
                           cost = cost,
                           method = "java",
                           parallel.java = cores.java, alloc_ram = ram.java)
    distance_result <- distance_result * res(resistance_1)[1]

    distance_result <- as.matrix(distance_result)
    rownames(distance_result)<- x[[id]]; colnames(distance_result)<- x[[id]]

    distance_result_matrix <- distance_result; distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
    distance_result <- as.data.frame(as.table(distance_result)); distance_result <- na.omit(distance_result)
    distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
    distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
    names(distance_result)<-c("From", "To", "Distance")

    #Inf
    x2 <- distance_result[which(is.infinite(distance_result[,3])),]

    if(nrow(x2) >= 1){
      Infinitos <- paste0(x2$From, "_", x2$To)
      distance_result2 <- distance_result[which(!is.infinite(distance_result[,3])),]

      error <- unique(x2$To); min_dist <- euclidean_distances(x, id, centroid = TRUE, distance_unit = "m")
      #
      if (is.null(geometry_out)){
        distance2 <- min_dist; distance2$idn <- paste0(distance2$From, "_", distance2$To)
        #Filter distances
        distance2 <- distance2[which(distance2$idn %in% Infinitos),]
        distance2$idn <- NULL; names(distance_result2)[3] <- "Distance"
        distance_result <- rbind(distance_result2, distance2)

        for(i in 1:nrow(distance2)){
          distance_result_matrix[which(rownames(distance_result_matrix) == distance2$From[i]),
                                 which(colnames(distance_result_matrix) == distance2$To[i])] <- distance2[i,3]

          distance_result_matrix[which(rownames(distance_result_matrix) == distance2$To[i]),
                                 which(colnames(distance_result_matrix) == distance2$From[i])] <- distance2[i,3]
        }

      } else {
        resistance_1 <- extend(resistance_1, as(x, 'Spatial')); resistance_1[is.na(resistance_1[])] <- geometry_out
        cost <- data.frame(code = geometry_out, cost = geometry_out)
        #Get new distances
        distance_result <- mat_cost_dist(raster = resistance_1,
                                         pts = pts_1,
                                         cost = cost,
                                         method = "java",
                                         parallel.java = cores.java, alloc_ram = ram.java)

        distance_result <- as.matrix(distance_result)
        rownames(distance_result)<- x[[id]]; colnames(distance_result)<- x[[id]]
        distance_result_matrix <- distance_result
        distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
        distance_result <- as.data.frame(as.table(distance_result)); distance_result <- na.omit(distance_result)
        distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
        distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
        names(distance_result)<-c("From", "To", "Distance")
      }
    }
  } else {
    if(is.null(bounding_circles)){
      st_geometry(coordenates_1) <- NULL
      coordenates_1 <- cbind(coordenates_1$lon, coordenates_1$lat)
      rownames(coordenates_1) <- NULL

      if(isTRUE(LCD)){
        Iso_conductance <- transition(resistance_1, CostFun, directions = ngh) %>%
          geoCorrection(., scl = FALSE)
        distance_result <- costDistance(Iso_conductance, coordenates_1) #reciprocal of conductance, works with resistance
      } else {
        Iso_conductance <- transition(resistance_1, CostFun, directions = ngh) %>%
          geoCorrection(., type = "r")
        distance_result <- commuteDistance(Iso_conductance, coordenates_1)
      }

      distance_result <- as.matrix(distance_result)
      rownames(distance_result)<- x[[id]]; colnames(distance_result)<- x[[id]]

      distance_result_matrix <- distance_result; distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
      distance_result <- as.data.frame(as.table(distance_result)); distance_result <- na.omit(distance_result)
      distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
      distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
      names(distance_result)<-c("From", "To", "Distance")

      #Inf
      x2 <- distance_result[which(is.infinite(distance_result[,3])),]

      if(nrow(x2) >= 1){
        Infinitos <- paste0(x2$From, "_", x2$To)
        distance_result2 <- distance_result[which(!is.infinite(distance_result[,3])),]

        error <- unique(x2$To); min_dist <- euclidean_distances(x, id, centroid = TRUE, distance_unit = "m")
        #
        if (is.null(geometry_out)){
          distance2 <- min_dist; distance2$idn <- paste0(distance2$From, "_", distance2$To)
          #Filter distances
          distance2 <- distance2[which(distance2$idn %in% Infinitos),]
          distance2$idn <- NULL; names(distance_result2)[3] <- "Distance"
          distance_result <- rbind(distance_result2, distance2)

          for(i in 1:nrow(distance2)){
            distance_result_matrix[which(rownames(distance_result_matrix) == distance2$From[i]),
                                   which(colnames(distance_result_matrix) == distance2$To[i])] <- distance2[i,3]

            distance_result_matrix[which(rownames(distance_result_matrix) == distance2$To[i]),
                                   which(colnames(distance_result_matrix) == distance2$From[i])] <- distance2[i,3]
          }

        } else {
          Iso_conductance <- raster(Iso_conductance) %>% extend(., as(x, 'Spatial'))

          Iso_conductance[is.na(Iso_conductance[])] <- geometry_out; Iso_conductance <- new("TransitionLayer", Iso_conductance)

          #Get new distances
          if(isTRUE(LCD)){
            distance_result<- costDistance(Iso_conductance, coordenates_1) #reciprocal of conductance, works with resistance
          } else {
            distance_result <- commuteDistance(Iso_conductance, coordenates_1)
          }

          distance_result <- as.matrix(distance_result)
          rownames(distance_result)<- x[[id]]; colnames(distance_result)<- x[[id]]

          distance_result_matrix <- distance_result

          distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
          distance_result <- as.data.frame(as.table(distance_result)); distance_result <- na.omit(distance_result)
          distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
          distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
          names(distance_result)<-c("From", "To", "Distance")
        }
      }
    } else {
      ####voronoi
      points_1 <- st_centroid(x); points_2 <- st_combine(st_geometry(points_1))
      voronoi_adj <- st_voronoi(points_2) %>% st_collection_extract(.)
      voronoi_adj <- voronoi_adj[unlist(st_intersects(points_1, voronoi_adj))] %>%
        st_touches(., sparse = FALSE) %>% as.matrix(.)
      colnames(voronoi_adj) <- x[[id]]; rownames(voronoi_adj) <- x[[id]]

      works <- as.numeric(availableCores())/2

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      plan(strategy = strat, gc = TRUE, workers = works)

      distance <- tryCatch(future_map(as.list(colnames(voronoi_adj)), function(i){
        foc_adj <- rbind(x[which(x[[id]] == i),],
                         x[as.vector(which(voronoi_adj[, which(x[[id]] == i)])),])
        foc_adj <- as(foc_adj, 'Spatial')
        dist <- euclidean_distances(foc_adj, id, threshold = bounding_circles, distance_unit = "m" )
        foc_adj <- foc_adj[which(foc_adj[[id]] %in% c(unique(dist$From), unique(dist$To))),]

        if(nrow(foc_adj) > 1){
          r2 <- gCentroid(foc_adj, byid = TRUE) %>% buffer(., width = bounding_circles) %>%
            crop(resistance_1, .)
          r2 <- raster::mask(r2, r2)

          ##
          coordenates_2 <- coordenates_1[which(coordenates_1[[id]] %in% foc_adj[[id]]),]
          st_geometry(coordenates_2) <- NULL
          coordenates_2 <- cbind(coordenates_2$lon, coordenates_2$lat)
          rownames(coordenates_2) <- NULL
          ##

          if(isTRUE(LCD)){
            Iso_conductance <- transition(r2, CostFun, directions = ngh) %>%
              geoCorrection(., scl = FALSE)
            distance_result <- costDistance(Iso_conductance, coordenates_2) #reciprocal of conductance, works with resistance
          } else {
            Iso_conductance <- transition(r2, CostFun, directions = ngh) %>%
              geoCorrection(., type = "r")
            distance_result <- commuteDistance(Iso_conductance, coordenates_2)
          }

          distance_result <- as.matrix(distance_result)
          rownames(distance_result)<- foc_adj[[id]]; colnames(distance_result)<- foc_adj[[id]]

          distance_result_matrix <- distance_result

          distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
          distance_result <- as.data.frame(as.table(distance_result));distance_result <- na.omit(distance_result)
          distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
          distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
          names(distance_result)<-c("From", "To", "Distance")

          #Inf
          x2 <- distance_result[which(is.infinite(distance_result[,3])),]

          if(nrow(x2) >= 1){
            Infinitos <- paste0(x2$From, "_", x2$To)
            distance_result2 <- distance_result[which(!is.infinite(distance_result[,3])),]

            error <- unique(x2$To)
            min_dist <- euclidean_distances(foc_adj, id, centroid = TRUE, distance_unit = "m")

            if (is.null(geometry_out)){
              distance2 <- min_dist;distance2$idn <- paste0(distance2$From, "_", distance2$To)
              #Filter distances
              distance2 <- distance2[which(distance2$idn %in% Infinitos),]
              distance2$idn <- NULL;names(distance_result2)[3] <- "Distance"
              distance_result <- rbind(distance_result2, distance2)

              for(i in 1:nrow(distance2)){
                distance_result_matrix[which(rownames(distance_result_matrix) == distance2$From[i]),
                                       which(colnames(distance_result_matrix) == distance2$To[i])] <- distance2[i,3]

                distance_result_matrix[which(rownames(distance_result_matrix) == distance2$To[i]),
                                       which(colnames(distance_result_matrix) == distance2$From[i])] <- distance2[i,3]
              }

            } else {
              Iso_conductance <- raster(Iso_conductance) %>%
                extend(., as(x, 'Spatial'))

              Iso_conductance[is.na(Iso_conductance[])] <- geometry_out; Iso_conductance <- new("TransitionLayer", Iso_conductance)

              #Get new distances
              if(isTRUE(LCD)){
                distance_result<- costDistance(Iso_conductance, coordenates_2) #reciprocal of conductance, works with resistance
              } else {
                distance_result <- commuteDistance(Iso_conductance, coordenates_2)
              }
              rm(Iso_conductance)

              distance_result <- as.matrix(distance_result)
              rownames(distance_result)<- foc_adj[[id]];colnames(distance_result)<- foc_adj[[id]]

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
          }
        }, .progress = TRUE), error = function(err) err)
      close_multiprocess(works)

      if (inherits(distance, "error")) {
        stop(distance)
      }

      distance <- purrr::compact(distance)

      if(isFALSE(pairwise)){
        uno <- lapply(1:length(distance), function(i){
          dos <- rownames(distance[[i]])
          return(dos)
        }) %>% do.call(c, .) %>% unique(.)

        mm <- matrix(nrow = length(uno), ncol = length(uno))
        rownames(mm) <- uno; colnames(mm) <- uno

        for(i in uno){
          mm2 <- mm[,i]; mm2 <- as.matrix(mm2)

          for(j in rownames(mm2)){
            if(j != i){
              mmj <- lapply(distance, function(x){
                x2 <- x[which(rownames(x)==i),which(rownames(x)==j)]
                if(length(x2) == 0){
                  x2 = NULL
                }
                return(x2)}) %>%
                purrr::compact(.) %>% unlist(.)
              mm2[j,] <- min(mmj)
            } else {
              mm2[j,] <- 0
            }
            mm[,i] <- mm2
          }
        }
        distance_result_matrix <- mm
      } else {
        distance_result <- do.call(rbind, distance); distance_result$dn <- paste0(distance_result$From, "_",distance_result$To)
        distance_result <- lapply(unique(distance_result$dn), function(i){
          d1 <- distance_result[which(distance_result$dn == i),];d1 <- d1[which(d1$Distance == min(d1$Distance)),1:3]
          if(nrow(d1)>1){
            d1 <- d1[1,]
          }
          return(d1)})
        distance_result <- do.call(rbind, distance_result);rownames(distance_result) <- NULL
      }
    }
  }

  if(isTRUE(pairwise)){
    if (!is.null(threshold)){
      distance_result <-  distance_result[which(distance_result[,3] <= threshold),]
    }
    rownames(distance_result) <- NULL
  } else {
    distance_result <- distance_result_matrix; distance_result <- as.matrix(distance_result)
    if (!is.null(threshold)){
      distance_result[which(distance_result <= threshold)] <- NA
    }
  }

  if(!is.null(write_table)){
    write.table(distance_result, write_table, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  return(distance_result)
}
