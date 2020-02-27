#' Get a Distance table that can be used as a connector file in Conefor
#'
#' Generates a distance table between nodes.
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. The shapefile must be in a projected coordinate system.
#' @param id character. Column name with the "node ID" in the shapefile data table.
#' @param type character. Choose one of the distances: "centroid" (faster, default), where Euclidean distance is
#' calculated from feature centroid; "edge", where Euclidean distance is calculated from feature edges;
#' "least-cost" that takes into account obstacles and local friction of the landscape (see, "gdistance" package);
#' "commute-time" that is analogous to the resistance distance of circuitscape. The commute-time distance is
#' based on the random walk theory and calculated using the electrical circuit theory (See, gdistance package).
#' If the type is equal to "least-cost" or "commute-time", then you have to use the "resistance" argument.
#' @param distance_unit character. If euclidean distance is selected you can set a distance unit, "Makurhini::unit_covert()"
#' compatible unit ("m", "km", "inch", "foot", "yard", "mile"). Default equal to meters "m".
#' @param keep numeric. Argument for higher processing speed. In case you have selected the "edge" distance, use this option to simplify the geometry and reduce the
#'  number of vertices (from rmapshaper::ms_simplif). The value can range from 0 to 1 and is the proportion of points to retain (default 0.02). The higher the value,
#'   the higher the speed but the greater uncertainty.
#' @param resistance raster. Raster object with resistance values.
#' @param CostFun A function to compute the cost to move between cells. Available only if you you have selected
#' the "least-cost" or "commute-time" distance. The default is the mean (isotropic cost distance):
#' function(x) mean(x).
#' @param ngh numeric.  Neighbour graph (directions) for distance calculations: 4 (von Neu-mann neighbourhood),
#' 8 (Moore neighbourhood) or 16. Available only if you you have selected
#' the "least-cost" or "commute-time" distance.
#' @param mask object of class sf, sfc, sfg or SpatialPolygons. For higher processing speed of "least-cost" or
#' "commute-time" distances. Use this option to clip the resistance at the extent of the mask.
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param geometry_out numeric. You can use this argument if you have selected the "least-cost" or "commute-time" distance.
#' If some spatial geometries are out of the resistance extent, then a buffer zone the large enough to cover
#' these spatial geometries and with this numeric value will be added to the resistance, so that it is possible to
#' calculate a cost distance value for the pairs of nodes that involve these geometries and and avoid an error.
#' If NULL, then a Euclidean distance (centroid) will be calculated to find these distances.
#' @param multiple character. Column name of the column with the zone to which each core area belongs.
#' @param prefix character. Initial prefix, use in case of processing several sites at the same time in CONEFOR command line.
#' @param write character. Choose the output folder if you use name_zona otherwise, place the output path, with
#' the name and extension ".txt".
#' @return Exports a euclidean or cost distance table between pairs of nodes.
#' @references
#' \url{https://www.rdocumentation.org/packages/rgeos/versions/0.3-26/topics/gDistance}\cr
#' \url{http://cgm.cs.mcgill.ca/~godfried/teaching/cg-projects/98/normand/main.html}.
#' @export
#' @importFrom methods as

distancefile <- function(nodes, id,
                         type =  "centroid",
                         distance_unit = NULL,
                         keep = NULL,
                         resistance = NULL,
                         CostFun = NULL,
                         ngh = NULL,
                         mask = NULL,
                         threshold = NULL,
                         geometry_out = NULL,
                         multiple = NULL, prefix = NULL,
                         write = NULL){

  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(nodes)[1] == "sf") {
    nodes <- as(nodes, 'Spatial')
    }

  if (is.null(distance_unit)) {
    distance_unit = "m"
    }

  if (nrow(nodes) > 1) {
    if (is.null(multiple)) {
      if (type %in%  c("centroid", "edge")){
        distance <- euclidean_distances(x = nodes, id = id, type_distance = type, distance_unit =distance_unit,
                                        keep = keep, threshold = threshold,  write_table = write)
        return(distance)

        } else if (type %in%  c("least-cost", "commute-time")){
          distance <- cost_distances(x = nodes, id = id, type_distance = type, resistance = resistance,
                                     CostFun = CostFun, ngh = ngh,
                                     threshold = threshold, mask = mask, geometry_out = geometry_out,
                                     write_table = write)
          return(distance)

          } else {
            stop("Error, you have to choose a type_distance option")
          }

    } else {
      multiple_1 <- unique(nodes@data[,which(colnames(nodes@data) == multiple)])
      if (type %in%  c("centroid", "edge")){
        distance <- lapply(as.list(multiple_1), function(x) { if (!is.null(write)){
          save <- paste(write, x, ".txt") } else { save <- NULL }
          nodes.1 <- nodes[nodes@data[,which(colnames(nodes@data) == multiple)] == x,]
          distance <- euclidean_distances(x = nodes.1, id = id, type_distance = type, distance_unit = distance_unit,
                                          keep = keep, threshold = threshold, write_table = save)
          return(distance) })
        names(distance) <- multiple_1
        return(distance)
        } else if (type %in%  c("least-cost", "commute-time")){
          distance <- lapply(as.list(multiple_1), function(x) { if (!is.null(write)){
            save <- paste(write, x, ".txt") } else { save <- NULL }
            nodes.1 <- nodes[nodes@data[,which(colnames(nodes@data) == multiple)] == x,]
            distance <- cost_distances(x = nodes.1, id = id, type_distance = type, resistance = resistance,
                                       CostFun = CostFun, ngh = ngh,
                                       threshold = threshold, mask = mask, geometry_out = geometry_out,
                                       write_table = save)
            return(distance) })
          names(distance) <- multiple_1
          return(distance)
          } else {
          stop("Error, you have to choose a type_distance option")
          }
      }
    } else {
      stop("Error, only one node")
    }
  }


