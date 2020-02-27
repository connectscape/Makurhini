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
#' @param write_table character. "" indicates output to the console.
#' @return Pairwise Euclidean distance table
#' @references Douglas, David H. and Peucker, Thomas K. (1973) "Algorithms for the Reduction of the Number of Points Required to Represent a Digitised Line or its Caricature", The Canadian Cartographer, 10(2), pp112-122.
#' @export
#' @importFrom rgeos gCentroid gDistance
#' @importFrom rmapshaper ms_simplify
#' @importFrom methods as
#' @importFrom utils combn write.table

euclidean_distances <- function(x, id, type_distance = "centroid", distance_unit = "m",
                               keep = NULL, threshold = NULL , write_table = NULL){
  if(missing(id)){
    stop("missing id")
  }
  if(missing(x)){
    stop("missing x object")
  }

  if(class(x)[1] == "sf") {
    x <- as(x, 'Spatial')
  }

  if (!is.null(write_table)) {
    if (!dir.exists(dirname(write_table))) {
      stop("error, output folder does not exist")
    }
  }

  if (!is.null(keep)){
    x_id <- x@data[,which(colnames(x@data) == id)]
    x <- ms_simplify(input = x, keep = keep, keep_shapes = TRUE, explode = TRUE)
    x$id <- x_id
    names(x) <- id
  }

  if (type_distance ==  "centroid"){
    centroid_1 <- gCentroid(x, byid = TRUE)
    distance <- gDistance(centroid_1, byid = TRUE)
    } else if (type_distance == "edge"){
      distance <- gDistance(x, byid = TRUE)
      } else if (type_distance == "hausdorff-edge"){
        distance <- gDistance(x, byid = TRUE, hausdorff = TRUE)
      } else {
        stop("Error, you have to choose a type_distance option")
      }

  name <- c(unique(x@data[,which(colnames(x@data) == id)]), "income")
  name <- name[1:(length(name)-1)]
  colnames(distance) <- name
  rownames(distance) <- name
  xy <- t(combn(colnames(distance), 2))
  distance_2 <- data.frame(xy, distance_2 = distance[xy])
  distance_2[,1] <- as.character(distance_2[,1])
  distance_2[,2] <- as.character(distance_2[,2])
  distance_2[,1] <- as.numeric(distance_2[,1])
  distance_2[,2] <- as.numeric(distance_2[,2])
  names(distance_2) <- c("From", "To", "Distance")

  if (distance_unit != "m"){
    distance_2$Distance <- unit_convert(data_unit = distance_2[,3], unit_1 = "m", unit_2 = distance_unit)
  }

  if (!is.null(threshold)){
    distance_2 <-  distance_2[which(distance_2$Distance <= threshold),]
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
#' @param CostFun A function to compute the cost to move between cells.The default is the mean (isotropic cost distance):  function(x) mean(x).
#' @param ngh numeric.  Neighbour graph (directions) for distance calculations: 4 (von Neu-mann neighbourhood), 8 (Moore neighbourhood) or 16 (king’s and knight’s moves). Default equal 16.
#' @param mask object of class sf, sfc, sfg or SpatialPolygons. For higher processing speed use this option to clip the resistance at the extent of the mask.
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param geometry_out numeric. If some spatial geometries are out of the resistance extent, then a buffer zone the large enough to cover these spatial geometries and with this numeric value will be added to the resistance, so that it is possibñe to calculate a cost distance value for the pairs of nodes that involve these geometries.
#'  If NULL, then a Euclidean distance will be calculated to find these distances.
#' @param write_table character.  "" indicates output to the console.
#' @return cost distance matrix
#' @references https://cran.r-project.org/web/packages/gdistance/gdistance.pdf
#' @export
#' @importFrom magrittr %>%
#' @import sf
#' @importFrom purrr map_dbl
#' @importFrom raster crop buffer maxValue mosaic
#' @importFrom gdistance transition geoCorrection costDistance commuteDistance
#' @importFrom methods as
#' @importFrom stats na.omit
#' @importFrom utils combn write.table

cost_distances <- function(x, id, type_distance = "least-cost", resistance = NULL, CostFun = NULL, ngh = NULL,
                           mask = NULL, threshold = NULL, geometry_out = NULL, write_table = NULL){

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

  st_geometry(coordenates_1)<- NULL
  coordenates_1 <- cbind(coordenates_1$lon, coordenates_1$lat)
  rownames(coordenates_1) <- NULL

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


  conductance <- transition(resistance_1, CostFun, directions = ngh)
  Iso_conductance <- geoCorrection(conductance, scl = FALSE)

  if(type_distance == "least-cost"){
    distance_result <- costDistance(Iso_conductance, coordenates_1) #reciprocal of conductance, works with resistance
    distance_result <- as.matrix(distance_result)
    rownames(distance_result)<- x[[id]]
    colnames(distance_result)<- x[[id]]
    distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
    distance_result <- as.data.frame(as.table(distance_result))
    distance_result <- na.omit(distance_result)
    distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
    distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
    names(distance_result)<-c("From", "To", "CostDistance")

  } else if (type_distance == "commute-time") {
    Iso_conductance <- geoCorrection(conductance, type = "r")
    distance_result <- commuteDistance(Iso_conductance, coordenates_1) #reciprocal of conductance, works with resistance
    distance_result <- as.matrix(distance_result)
    rownames(distance_result)<- x[[id]]
    colnames(distance_result)<- x[[id]]
    distance_result[lower.tri(distance_result, diag = TRUE)] <- NA
    distance_result <- as.data.frame(as.table(distance_result))
    distance_result <- na.omit(distance_result)
    distance_result[,1] <- as.numeric(as.character(distance_result[,1]))
    distance_result[,2]<-as.numeric(as.character(distance_result[,2]))
    names(distance_result)<-c("From", "To", "CommuteDistance")
    } else {
      stop("Error, you have to choose a type_distance option")
    }

  x2 <- distance_result[which(is.infinite(distance_result[,3])),]#infinitos
  Infinitos <- paste0(x2$From, "_", x2$To)
  distance_result2 <- distance_result[which(!is.infinite(distance_result[,3])),]

  if(nrow(x2) >= 1){
    error <- unique(x2$To)
    unique(x2$From)
    x2 <- x[error,]
    x4 <- x
    min_dist <- euclidean_distances(x4, id, type_distance = "centroid", distance_unit = "m")
    #
    if (is.null(geometry_out)){
      distance2 <- min_dist
      distance2$idn <- paste0(distance2$From, "_", distance2$To)
      #Filter distances
      distance2 <- distance2[which(distance2$idn %in% Infinitos),]
      distance2$idn <- NULL
      names(distance_result2)[3] <- "Distance"
      distance_result <- rbind(distance_result2, distance2)
      } else {
       #Encontrar distancia buffer
        minimum <- min_dist
        minimum$idn <- paste0(minimum$From, "_", minimum$To)
        minimum <- minimum[which(minimum$idn %in% Infinitos),]
        minimum$idn <- NULL
        minimum <- max(minimum)
        #Buffer
       rbuffer <- buffer(resistance_1, width = minimum, doEdge = TRUE)
       resistance_2 <- resistance_1
       resistance_2[resistance_2 > 0] <- maxValue(resistance_1) + 10
       resistance_2 <- mosaic(resistance_2, rbuffer, fun = max)
       resistance_2[resistance_2 == 1] <- maxValue(resistance_1) + 20
       #Mosaic
       resistance_2 <- mosaic(resistance_2, resistance_1, fun = min)
       resistance_2[resistance_2 ==  maxValue(resistance_1) + 20] <- geometry_out

       #Conductancia
       conductance <- transition(resistance_2, CostFun, directions = ngh)
       Iso_conductance <- geoCorrection(conductance)

       coordenates_2 <- x4  %>%
         mutate(lon = map_dbl(geometry, ~st_centroid_within_poly(.x)[[1]]),
                       lat = map_dbl(geometry, ~st_centroid_within_poly(.x)[[2]]))
       st_geometry(coordenates_2)<- NULL
       coordenates_2 <- cbind(coordenates_2$lon, coordenates_2$lat)
       rownames(coordenates_2) <- NULL
       coordenates_2 <- rbind(coordenates_2, coordenates_1)

       #Get new distances
       if(type_distance == "least-cost"){
         distance2 <- costDistance(Iso_conductance, coordenates_2) #reciprocal of conductance, works with resistance
         distance2 <- as.matrix(distance2)
         rownames(distance2)<- x4[[id]]
         colnames(distance2)<- x4[[id]]
         distance2[lower.tri(distance2, diag = TRUE)] <- NA
         distance2 <- as.data.frame(as.table(distance2))
         distance2 <- na.omit(distance2)
         distance2[,1] <- as.numeric(as.character(distance2[,1]))
         distance2[,2]<-as.numeric(as.character(distance2[,2]))
         names(distance2)<-c("From", "To", "CostDistance")
       } else if (type_distance == "commute-time") {
         Iso_conductance <- geoCorrection(conductance, type = "r")
         distance2 <- commuteDistance(Iso_conductance, coordenates_2) #reciprocal of conductance, works with resistance
         distance2 <- as.matrix(distance2)
         rownames(distance2)<- x4[[id]]
         colnames(distance2)<- x4[[id]]
         distance2[lower.tri(distance2, diag = TRUE)] <- NA
         distance2 <- as.data.frame(as.table(distance2))
         distance2 <- na.omit(distance2)
         distance2[,1] <- as.numeric(as.character(distance2[,1]))
         distance2[,2]<-as.numeric(as.character(distance2[,2]))
         names(distance2)<-c("From", "To", "CommuteDistance")
       }
       #Filter distances
       distance2$idn <- paste0(distance2$From, "_", distance2$To)
       #Filter distances
       distance2 <- distance2[which(distance2$idn %in% Infinitos),]
       distance2$idn <- NULL

       a <- paste0(distance2$From, "_", distance2$To)
       distance2 <- distance2[which(!duplicated(a)),]

       #Union
       distance_result <- rbind(distance_result2, distance2)
      }
  }


  if (!is.null(threshold)){
    distance_result <-  distance_result[which(distance_result[,3] <= threshold),]
  }

  if(!is.null(write_table)){
    write.table(distance_result, write_table, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  return(distance_result)
  }
