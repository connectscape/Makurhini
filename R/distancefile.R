#' Get a Distance table that can be used as a connector file in Conefor
#'
#' Generates a distance table between nodes.
#' @param nodes object of class sf, SpatialPolygons, RasterLaryer or SpatRaster (terra package). It must be in a projected coordinate system.
#' If nodes is a raster layer then raster values (Integer) will be taken as "id".
#' @param id character. If nodes is a shapefile then you must specify the column name with the node ID in the shapefile data table.
#' @param type character. Choose one of the distances: "centroid" (faster, default), where Euclidean distance is
#' calculated from feature centroid; "edge", where Euclidean distance is calculated from feature edges;
#' "least-cost" that takes into account obstacles and local friction of the landscape (see, "gdistance" package);
#' "commute-time" that is analogous to the resistance distance of circuitscape. The commute-time distance is
#' based on the random walk theory and calculated using the electrical circuit theory (See, gdistance package).
#' If the type is equal to "least-cost" or "commute-time", then you have to use the "resistance" argument.
#' @param distance_unit character. If Euclidean distance is selected you can set a distance unit, "Makurhini::unit_covert()"
#' compatible unit ("m", "km", "inch", "foot", "yard", "mile"). Default equal to meters "m".
#' @param keep numeric. Argument for higher processing speed. In case you have selected the "edge" distance, use this option to simplify the geometry and reduce the
#'  number of vertices (from rmapshaper::ms_simplify). The value can range from 0 to 1 and is the proportion of points to retain (default equal to NULL). The higher the value,
#'   the higher the speed but less precision.
#' @param resistance raster. Raster object with resistance values.
#' @param CostFun A function to compute the cost to move between cells. Available only if you you have selected
#' the "least-cost" or "commute-time" distance. The default is the mean (isotropic cost distance):
#' function(x) mean(x). If resistance is a conductance raster then you can use: function(x) 1/mean(x)
#' @param ngh numeric.  Neighbor graph (directions) for distance calculations: 4 (von Neu-mann neighbourhood),
#' 8 (Moore neighborhood) or 16. Available only if you have selected
#' the "least-cost" or "commute-time" distance.
#' @param mask object of class sf, sfc, sfg or SpatialPolygons. For higher processing speed of "least-cost" or
#' "commute-time" distances. Use this option to clip the resistance at the extent of the mask.
#' @param threshold numeric. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param geometry_out numeric. You can use this argument if you have selected the "least-cost" or "commute-time" distance.
#' If some spatial geometries are out of the resistance extent, then a buffer zone the large enough to cover
#' these spatial geometries and with this numeric value will be added to the resistance, so that it is possible to
#' calculate a cost distance value for the pairs of nodes that involve these geometries and avoid an error.
#' If NULL, then a Euclidean distance (centroid) will be calculated to find these distances.
#' @param multiple character. If nodes is a shapefile then you can use this argument. Name of the column with the zone to which each core area belongs.
#' @param prefix character. Initial prefix, use in case of processing several sites at the same time in CONEFOR command line.
#' @param parallel logical. If nodes is a raster then you can use this argument for larges RasterLayer.
#' @param bounding_circles numeric. If a value is entered, this will create bounding circles around pairs of core areas
#'  (recommended for speed, large resistance rasters or pixel resolution < 150 m).
#'  Buffer distances are entered in map units. Also, the function is parallelized using
#'   and furrr package and multiprocess plan, default = NULL.
#' @param pairwise logical. If TRUE a pairwise table is returned (From, To, distance) otherwise it will be a matrix.
#' @param write character. Choose the output folder if you use multiple argument otherwise, place the output path, with
#' the name and extension ".txt".
#' @return Exports a euclidean or cost distance table between pairs of nodes.
#' @references
#' \url{https://www.rdocumentation.org/packages/rgeos/versions/0.3-26/topics/gDistance}\cr
#' \url{http://cgm.cs.mcgill.ca/~godfried/teaching/cg-projects/98/normand/main.html}.
#' @export
#' @import sf
#' @importFrom methods as
#' @importFrom purrr map
#' @importFrom future multiprocess plan
#' @importFrom furrr future_map
#' @importFrom rmapshaper ms_dissolve
#' @importFrom spex qm_rasterToPolygons
#' @importFrom raster rasterToPoints crs values raster
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize

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
                         bounding_circles = NULL,
                         multiple = NULL, prefix = NULL,
                         parallel = FALSE, edgeParallel = FALSE,
                         pairwise = TRUE,
                         write = NULL){

  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if (is.null(distance_unit)) {
    distance_unit = "m"
  }

  if(class(nodes)[1] == "sf" | class(nodes)[1] == "SpatialPolygonsDataFrame"){

  if(class(nodes)[1] == "sf") {
    nodes <- as(nodes, 'Spatial')
    }

  if (nrow(nodes) > 1) {
    if (is.null(multiple)) {
      if (type %in%  c("centroid", "edge")){
        distance <- euclidean_distances(x = nodes, id = id, type_distance = type, distance_unit =distance_unit,
                                        keep = keep, threshold = threshold,  pairwise = pairwise,
                                        write_table = write)
        return(distance)

        } else if (type %in%  c("least-cost", "commute-time")){
          distance <- cost_distances(x = nodes, id = id, type_distance = type, resistance = resistance,
                                     CostFun = CostFun, ngh = ngh, bounding_circles = bounding_circles ,
                                     threshold = threshold, mask = mask, geometry_out = geometry_out,
                                     pairwise = pairwise,
                                     write_table = write)
          return(distance)

          } else {
            stop("Error, you have to choose a type_distance option")
          }

    } else {
      multiple_1 <- unique(nodes@data[,which(colnames(nodes@data) == multiple)])
      if (type %in%  c("centroid", "edge")){
        distance <- map(as.list(multiple_1), function(x) { if (!is.null(write)){
          save <- paste(write, x, ".txt") } else { save <- NULL }
          nodes.1 <- nodes[nodes@data[,which(colnames(nodes@data) == multiple)] == x,]
          distance <- euclidean_distances(x = nodes.1, id = id, type_distance = type, distance_unit = distance_unit,
                                          keep = keep, threshold = threshold, pairwise = pairwise,
                                          write_table = save)
          return(distance) })
        names(distance) <- multiple_1
        return(distance)
        } else if (type %in%  c("least-cost", "commute-time")){
          distance <- map(as.list(multiple_1), function(x) { if (!is.null(write)){
            save <- paste(write, x, ".txt") } else { save <- NULL }
            nodes.1 <- nodes[nodes@data[,which(colnames(nodes@data) == multiple)] == x,]
            distance <- cost_distances(x = nodes.1, id = id, type_distance = type, resistance = resistance,
                                       CostFun = CostFun, ngh = ngh, bounding_circles = bounding_circles,
                                       threshold = threshold, mask = mask, geometry_out = geometry_out,
                                       pairwise = pairwise,
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

  } else if(class(nodes)[1] == "RasterLayer"| class(nodes)[1] == "SpatRaster"){
    if(class(nodes)[1] == "SpatRaster"){
      nodes <- raster(nodes)
    }

    if (nodes@data@max > 1) {
      if (is.null(multiple)) {
        if (type %in%  c("centroid", "edge")){
          if(type == "centroid"){
            if(isTRUE(parallel)){
                rp <- unique(raster::values(nodes))
                rp <- as.vector(rp)
                rp <- rp[which(!is.na(rp))]
                rp <- split(rp, ceiling(seq_along(rp)/round(sqrt(max(rp)))))

                '%!in%' <- function(x,y)!('%in%'(x,y))
                plan(strategy = multiprocess)

              coords <- future_map(rp, function(x){
                r <- nodes
                r[r %!in% x] <- NA
                p <- data.frame(rasterToPoints(r))
                names(p)[3] <- "layer"
                p <- p[p$layer > 0,]

                coords1 <- map(split(p, p$layer), function(x){
                  cn <- colMeans(x[, c("x", "y")])
                  cn <- data.frame(x = cn[[1]], y = cn[[2]], Id = unique(x[, "layer"]))
                  return(cn)}) %>% do.call(rbind, .)

                coords1 <- st_as_sf(coords1, coords = c("x", "y"),
                                   crs = crs(nodes), stringsAsFactors = FALSE)
                return(coords1) },.progress = TRUE)
              coords <- do.call(rbind, coords)


            } else {
              p <- data.frame(rasterToPoints(nodes))
              names(p)[3] <- "layer"

              p <- p[p$layer > 0,]

              coords <- map(split(p, p$layer), function(x){
                cn <- colMeans(x[, c("x", "y")])
                cn <- data.frame(x = cn[[1]], y = cn[[2]], Id = unique(x[, "layer"]))
                return(cn)}) %>% do.call(rbind, .)

              coords <- st_as_sf(coords, coords = c("x", "y"),
                                 crs = crs(nodes), stringsAsFactors = FALSE)
            }

            distance <- euclidean_distances(x = coords, id = "Id", type_distance = "centroid",
                                            distance_unit = distance_unit,
                                            keep = keep, threshold = threshold,
                                            pairwise = pairwise,
                                            write_table = write)

          } else {
            if(isTRUE(parallel)){
              rp <- unique(values(nodes))
              rp <- as.vector(rp)
              rp <- rp[which(!is.na(rp))]
              rp <- split(rp, ceiling(seq_along(rp)/round(sqrt(max(rp)))))

              '%!in%' <- function(x,y)!('%in%'(x,y))
              plan(strategy = multiprocess)

              pol_nodes <- future_map(rp, function(x){
                r <- nodes
                r[r %!in% x] <- NA
                pn <- qm_rasterToPolygons(r)
                names(pn)[1] <- "Id"
                pn <- pn %>% group_by(.data$Id) %>% summarize() %>% st_cast(.)
                return(pn)
                }, .progress = TRUE) %>% do.call(rbind, .)



            } else {
              pol_nodes <- qm_rasterToPolygons(nodes)
              names(pol_nodes)[1] <- "Id"
              pol_nodes <- map(split(pol_nodes, pol_nodes$Id), function(x){
                x2 <- st_union(x) %>% st_sf(.)
                x2$Id <- unique(x$Id)
                return(x2)}) %>% do.call(rbind, .)

            }
            distance <- euclidean_distances(x = pol_nodes, id = "Id", type_distance = "edge",
                                            distance_unit = distance_unit, edgeParallel = parallel,
                                            keep = keep, threshold = threshold, pairwise = pairwise,

                                            write_table = write)

          }

          return(distance)

        } else if (type %in%  c("least-cost", "commute-time")){

          if(isTRUE(parallel)){
            rp <- unique(values(nodes))
            rp <- as.vector(rp)
            rp <- rp[which(!is.na(rp))]
            rp <- split(rp, ceiling(seq_along(rp)/round(sqrt(max(rp)))))

            '%!in%' <- function(x,y)!('%in%'(x,y))
            plan(strategy = multiprocess)

            coords <- future_map(rp, function(x){
              r <- nodes
              r[r %!in% x] <- NA
              p <- data.frame(rasterToPoints(r))
              p <- p[p$layer > 0,]

              coords1 <- map(split(p, p$layer), function(x){
                cn <- colMeans(x[, c("x", "y")])
                cn <- data.frame(x = cn[[1]], y = cn[[2]], Id = unique(x[, "layer"]))
                return(cn)}) %>% do.call(rbind, .)

              coords1 <- st_as_sf(coords1, coords = c("x", "y"),
                                  crs = crs(nodes), stringsAsFactors = FALSE)
              return(coords1) },.progress = TRUE)
            coords <- do.call(rbind, coords)


          } else {
            p <- data.frame(rasterToPoints(nodes))
            names(p)[3] <- "layer"

            p <- p[p$layer > 0,]

            coords <- map(split(p, p$layer), function(x){
              cn <- colMeans(x[, c("x", "y")])
              cn <- data.frame(x = cn[[1]], y = cn[[2]], Id = unique(x[, "layer"]))
              return(cn)}) %>% do.call(rbind, .)

            coords <- st_as_sf(coords, coords = c("x", "y"),
                               crs = crs(nodes), stringsAsFactors = FALSE)
          }

          distance <- cost_distances(x = coords, id = "Id", type_distance = type, resistance = resistance,
                                     CostFun = CostFun, ngh = ngh, bounding_circles = bounding_circles ,
                                     threshold = threshold, mask = mask, geometry_out = geometry_out,
                                     pairwise = pairwise,
                                     write_table = write)
          return(distance)

        } else {
          stop("Error, you have to choose a type_distance option")
        }

      }

    } else {
      stop("Error, only one node")
    }
  } else {
    stop("error missing file of nodes or check the class")
  }

  }


