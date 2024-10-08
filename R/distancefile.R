#' Get a table or matrix with the distances between pairs of nodes.
#'
#' Get a table or matrix with the distances (Euclidean or cost distances) between pairs of nodes. It can be used as a connector file.
#' @param nodes object of class \code{sf, SpatialPolygons, RasterLaryer or SpatRaster} (terra package). It must be in a \bold{projected coordinate system}.
#' If nodes is a \code{raster} layer then raster values (Integer) will be taken as "id".
#' @param id \code{character}. If nodes is a shapefile then you must specify the column name with the node ID in the shapefile data table.
#' @param type \code{character}. Choose one of the distance options:\cr
#' -    \code{"centroid"} (faster option, default), where Euclidean distance is calculated from feature centroid
#' -    \code{"edge"}, where Euclidean distance is calculated from feature edges.\cr
#' -    \code{"least-cost"} that takes into account obstacles and local friction of the landscape (see, \bold{gdistance} package).\cr
#' -    \code{"commute-time"} that is analogous to the resistance distance of circuitscape. The commute-time distance is
#' based on the random walk theory and calculated using the electrical circuit theory (See, \bold{gdistance} package).\cr
#' If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you have to use the \bold{"resistance"} argument.
#' @param distance_unit \code{character}. If \code{"least-cost"} or \code{"commute-time"} is selected or \code{resist.units = TRUE} then you can set a distance unit
#'  ("m", "km", "inch", "foot", "yard", "mile"). Default equal to meters "m" (see, \link[Makurhini]{unit_convert}).
#' @param keep \code{numeric}. In case you have selected the \code{"edge"} distance, use this option to simplify the geometry and reduce the
#'  number of vertices. The value can range from 0 to 1 and is the proportion of points to retain (default equal to \code{NULL}). The higher the value,
#'   the higher the speed but less precision.
#' @param resistance \code{raster}. Raster object with \bold{resistance} values. If \bold{\code{least_cost.java = TRUE}}, then
#'  resistance must bee an \bold{integer raster} (i.e., integer values).
#' @param resist.units \code{logical}. If \code{resist.units = TRUE} and \code{type = "least-cost"} then cost units are converted to metric units by multiplying the cost by the raster resolution.
#' @param CostFun \code{function}. A function to compute the cost to move between cells. Available only if you you have selected
#' the \code{"least-cost"} or \code{"commute-time"} distance. The default is the mean (isotropic cost distance):
#' \code{function(x) 1/mean(x)}. If resistance is a conductance raster then you can use: \code{function(x) mean(x)}
#' @param ngh \code{numeric}. Neighbor graph (directions) for distance calculations: 4 (von Neu-mann neighbourhood),
#' 8 (Moore neighborhood) or 16. Available only if you have selected the \code{"least-cost"} or \code{"commute-time"} distance.
#' @param mask object of class \code{sf, sfc, sfg. SpatialPolygons}. For higher processing speed of \code{"least-cost"} or \code{"commute-time"} distances. Use this option to clip the resistance at the extent of the mask.
#' @param threshold \code{numeric}. Distance threshold, pairs of nodes with a distance value above this threshold will be discarded.
#' @param geometry_out \code{numeric}. You can use this argument if you have selected the \code{"least-cost"} or \code{"commute-time"} distance.
#' If some spatial geometries are out of the resistance extent, then a buffer zone the large enough to cover
#' these spatial geometries and with this numeric value will be added to the resistance, so that it is possible to
#' calculate a cost distance value for the pairs of nodes that involve these geometries and avoid an error.
#' If NULL, then the Euclidean distance between centroids (\code{type = "centroid"}) will be calculated to find these distances.
#' @param parallel \code{logical} or \code{numerical}. Recommended for a large number of nodes or very large RasterLayer.
#' @param ActiveParallel \code{logical}. It should be \code{TRUE} if there is already an open parallelization plan.
#' @param bounding_circles \code{numeric}. If a value is entered, this will create bounding circles around pairs of core areas
#'  (recommended for speed, large resistance rasters or pixel resolution < 150 m).
#'  Buffer distances are entered in map units.
#' @param least_cost.java \code{logical}. If \code{TRUE} then the programming language and computing platform \bold{'java'} will be used
#' to estimate the least-cost distance \bold{using a resistance raster and the following formula}: \code{function(x) 1/mean(x)}. It is necessary to have \bold{java} installed. This option use
#' the package \bold{\code{graph4lg}} to reduce computation times.
#' @param cores.java \code{numeric}. Computer cores used to run the .jar file (see, \bold{\code{graph4lg}} ), default = \code{1}.
#' @param ram.java \code{numeric}. RAM gigabytes to run the .jar file (see, \bold{\code{graph4lg}} ), default = \code{NULL}.
#' @param pairwise \code{logical}. If \code{TRUE} a pairwise table of class \code{data.frame} is returned (From, To, distance) otherwise it will be a \code{matrix}.
#' @param write \code{character}. Output path, with name and extension \bold{".txt"}.
#' @return Exports a euclidean or cost distance table between pairs of nodes.
#' @references
#' \url{https://www.rdocumentation.org/packages/rgeos/versions/0.3-26/topics/gDistance}\cr
#' \url{http://cgm.cs.mcgill.ca/~godfried/teaching/cg-projects/98/normand/main.html}.
#' @export
#' @importFrom sf st_as_sf st_union st_sf
#' @importFrom methods as
#' @importFrom purrr map_dfr
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map_dfr
#' @importFrom magrittr %>%
#' @importFrom terra rast minmax as.polygons unique subst
#' @importFrom raster rasterToPoints crs raster
distancefile <- function(nodes, id, type =  "centroid", distance_unit = NULL,
                         keep = NULL,
                         resistance = NULL, resist.units = FALSE,
                         CostFun = NULL, ngh = NULL, mask = NULL,
                         threshold = NULL, geometry_out = NULL, bounding_circles = NULL,
                         parallel = FALSE, ActiveParallel = FALSE,
                         least_cost.java = FALSE,
                         cores.java = 1, ram.java = NULL,
                         pairwise = TRUE,
                         write = NULL){
  . = NULL
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

  if(!grepl("Raster", class(nodes)[1])){
    if(class(nodes)[1] != "sf" ){
      nodes <- st_as_sf(nodes)
    }

    if(nrow(nodes) < 2){
      stop("error, you need more than 2 nodes")
    }

    '%!in%' <- function(x,y)!('%in%'(x,y))

    if (is.null(id)) {
      stop("error missing id name")
    } else {
      if(length(which(names(nodes) == id)) == 0){
        stop("error missing id name in nodes file")
      }
    }

    if(isTRUE(ActiveParallel)){
      parallel = FALSE
      if(type == "edge"){
        message("Parallel will not be used to get the edges of the fragments but it will be used to get the distances that is because we depend on the 'terra' package which until the release of this version of Makurhini cannot be serialized. A sequential loop will be used.")
      }
    } else {
      if(isTRUE(parallel) | !is.null(parallel)){
        parallel = TRUE; works <- as.numeric(availableCores())-1

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        plan(strategy = strat, gc = TRUE, workers = works)
      }
    }
  } else {
    ###POR AHORA
    if(class(nodes)[1] != "RasterLayer"){
      nodes <- raster::raster(nodes)
    }

    if(isTRUE(ActiveParallel)){
      parallel = FALSE
      if(type == "edge"){
        message("Parallel will not be used to get the edges of the fragments but it will be used to get the distances that is because we depend on the 'terra' package which until the release of this version of Makurhini cannot be serialized. A sequential loop will be used.")
      }
    } else {
      if(isTRUE(parallel) | !is.null(parallel)){
        if(type != "edge"){
          if(is.numeric(parallel)){
            works <- parallel
          } else {
            works <- as.numeric(availableCores())-1
          }

          if(.Platform$OS.type == "unix") {
            strat <- future::multicore
          } else {
            strat <- future::multisession
          }
          plan(strategy = strat, gc = TRUE, workers = works)
        } else {
          message("Parallel will not be used to get the edges of the fragments but it will be used to get the distances that is because we depend on the 'terra' package which until the release of this version of Makurhini cannot be serialized. A sequential loop will be used.")
        }
      }
    }

    if(type == "centroid"){
      if(class(nodes)[1] == "SpatRaster"){
        nodes <- raster(nodes)
      }

      if(nodes@data@max < 2) {
        stop("error, you need more than 2 nodes")
      }

      if(isTRUE(parallel) | !is.null(parallel)){
        cr <- crs(nodes); pts <- data.frame(rasterToPoints(nodes)); pts <- pts[pts[[3]] > 0,]
        nodes <- tryCatch(future_map_dfr(split(pts, pts[[3]]), function(x){
          x.1 <- colMeans(x[, c("x", "y")]); x.1 <- data.frame("x" = x.1[[1]], "y" = x.1[[2]], "Id" = unique(x[[3]]))
          return(x.1)}), error = function(err)err)

        if(isFALSE(ActiveParallel)){
          close_multiprocess()
        }

        if(inherits(nodes, "error")) {
          close_multiprocess(); stop(nodes)
        } else {
          nodes <- st_as_sf(nodes, coords = c("x", "y"),  crs = cr, stringsAsFactors = FALSE)
        }
      } else {
        cr <- crs(nodes); pts <- data.frame(rasterToPoints(nodes)); pts <- pts[pts[[3]] > 0,]
        nodes <- map_dfr(split(pts, pts[[3]]), function(x){
          x.1 <- colMeans(x[, c("x", "y")]); x.1 <- data.frame("x" = x.1[[1]], "y" = x.1[[2]], "Id" = unique(x[[3]]))
          return(x.1)}, .progress = TRUE)
        nodes <- st_as_sf(nodes, coords = c("x", "y"),  crs = cr, stringsAsFactors = FALSE)
      }
    }

    if(type == "edge"){
      if(class(nodes)[1] != "SpatRaster"){
        nodes <- rast(nodes)
      }

      if(minmax(nodes)[2] < 2) {
        stop("error, you need more than 2 nodes")
      }

      if(isTRUE(parallel) | !is.null(parallel)){
        rp <- terra::unique(nodes, na.rm = TRUE)[[1]]; rp2 <- base::split(rp, ceiling(seq_along(rp)/round(sqrt(max(rp)))))

        nodes <- tryCatch(map_dfr(rp2, function(x){
          x.1 <- subst(nodes, rp[!rp %in% x], NA)
          x.1 <- suppressWarnings(as.polygons(x.1, na.all = FALSE) %>% st_as_sf(.))
          return(x.1)}, .progress = TRUE), error = function(err)err)

        if(inherits(nodes, "error")) {
          stop(nodes)
        }
      } else {
        nodes <- terra::as.polygons(nodes) %>% st_as_sf(.)
      }
    }
    names(nodes)[1] <- "IdTemp"; id = "IdTemp"
  }

  if(type %in%  c("centroid", "edge")){
    distance <- euclidean_distances(x = nodes,
                                    id = id,
                                    centroid = if(type == "centroid"){TRUE} else {FALSE},
                                    distance_unit = distance_unit,
                                    keep = keep, threshold = threshold,
                                    distParallel = parallel,
                                    ActiveParallel = ActiveParallel,
                                    pairwise = pairwise,
                                    write_table = write)
    return(distance)
  } else if (type %in%  c("least-cost", "commute-time")){
    distance <- cost_distances(x = nodes,
                               id = id,
                               LCD = if(type == "least-cost"){TRUE} else {FALSE},
                               resistance = resistance,
                               CostFun = CostFun,
                               ngh = ngh,
                               bounding_circles = bounding_circles ,
                               threshold = threshold,
                               mask = mask,
                               geometry_out = geometry_out,
                               pairwise = pairwise,
                               least_cost.java = least_cost.java,
                               cores.java = cores.java, ram.java = ram.java,
                               write_table = write)
    if(isTRUE(resist.units) & type == "least-cost"){
      distance <- distance * res(resistance)
      if(!is.null(distance_unit)){
        distance <- unit_convert(distance, "m", distance_unit)
      }
    }
    invisible(gc()); return(distance)
  } else {
    stop("Error, you have to choose a type option")
  }

}


