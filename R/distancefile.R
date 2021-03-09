#' Get a table or matrix with the distances between pairs of nodes.
#'
#' Get a table or matrix with the distances (Euclidean or cost distances) between pairs of nodes. It can be used as a connector file.
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
#'  number of vertices. The value can range from 0 to 1 and is the proportion of points to retain (default equal to NULL). The higher the value,
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
#' @param parallel logical. If nodes is a raster then you can use this argument for larges RasterLayer.
#' @param edgeParallel logical. Parallelize the edge distance using furrr package and multiprocess plan, default = FALSE.
#' @param bounding_circles numeric. If a value is entered, this will create bounding circles around pairs of core areas
#'  (recommended for speed, large resistance rasters or pixel resolution < 150 m).
#'  Buffer distances are entered in map units. Also, the function is parallelized using
#'   and furrr package and multiprocess plan, default = NULL.
#' @param least_cost.java logical. If TRUE then the programming language and computing platform 'java' will be used
#' to estimate the least-cost distance USING only the FORMULA: function(x) 1/mean(x). It is necessary to have java installed. This option use
#' the package 'graph4lg' to reduce computation times.
#' @param cores.java numeric. Computer cores used to run the .jar file (see, graph4lg), default = 1.
#' @param ram.java numeric. RAM gigabytes to run the .jar file (see, graph4lg), default = NULL.
#' @param pairwise logical. If TRUE a pairwise table is returned (From, To, distance) otherwise it will be a matrix.
#' @param write character. Output path, with name and extension ".txt".
#' @return Exports a euclidean or cost distance table between pairs of nodes.
#' @references
#' \url{https://www.rdocumentation.org/packages/rgeos/versions/0.3-26/topics/gDistance}\cr
#' \url{http://cgm.cs.mcgill.ca/~godfried/teaching/cg-projects/98/normand/main.html}.
#' @export
#' @importFrom sf st_as_sf st_cast st_union st_sf
#' @importFrom methods as
#' @importFrom purrr map
#' @importFrom future multiprocess plan availableCores
#' @importFrom furrr future_map
#' @importFrom rmapshaper ms_dissolve
#' @importFrom spex qm_rasterToPolygons
#' @importFrom raster rasterToPoints crs values raster
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize
#' @importFrom rlang .data
distancefile <- function(nodes, id, type =  "centroid", distance_unit = NULL, keep = NULL,
                         resistance = NULL, CostFun = NULL, ngh = NULL, mask = NULL,
                         threshold = NULL, geometry_out = NULL, bounding_circles = NULL,
                         parallel = FALSE, edgeParallel = FALSE,
                         least_cost.java = FALSE,
                         cores.java = 1, ram.java = NULL,
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

  if(class(nodes)[1] != "RasterLayer"){
    if(class(nodes)[1] == "sf" ){
      nodes <- as(nodes, 'Spatial')
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

  } else {
    if(class(nodes)[1] == "SpatRaster"){
      nodes <- raster(nodes)
    }
    if(nodes@data@max < 2) {
      stop("error, you need more than 2 nodes")
    }
  }

  #############################
  if(class(nodes)[1] == "RasterLayer"){
    if(isTRUE(parallel)){
      rp <- unique(values(nodes))
      rp <- as.vector(rp)
      rp <- rp[which(!is.na(rp))]
      rp <- split(rp, ceiling(seq_along(rp)/round(sqrt(max(rp)))))
      works <- as.numeric(availableCores())-1
      plan(strategy = multiprocess, gc = TRUE, workers = works)

      if(type == "centroid"){
        cr <- crs(nodes)
        nodes <- future_map(rp, function(x){
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
                              crs = cr, stringsAsFactors = FALSE)
          return(coords1) },.progress = TRUE) %>% do.call(rbind, .)
      } else {
        nodes <- future_map(rp, function(x){
          r <- nodes
          r[r %!in% x] <- NA
          pn <- qm_rasterToPolygons(r)
          names(pn)[1] <- "Id"
          pn <- pn %>% group_by(.data$Id) %>% summarize() %>% st_cast(.)
          return(pn)
        }, .progress = TRUE) %>% do.call(rbind, .)
      }
      close_multiprocess(works)

    } else {
      if(type == "centroid"){
        cr <- crs(nodes)
        pts <- data.frame(rasterToPoints(nodes))
        names(pts)[3] <- "layer"

        pts <- pts[pts$layer > 0,]

        nodes <- lapply(split(pts, pts$layer), function(x){
          cn <- colMeans(x[, c("x", "y")])
          cn <- data.frame(x = cn[[1]], y = cn[[2]], Id = unique(x[, "layer"]))
          return(cn)}) %>% do.call(rbind, .)

        nodes <- st_as_sf(nodes, coords = c("x", "y"),
                          crs = cr, stringsAsFactors = FALSE)
      } else {
        nodes <- qm_rasterToPolygons(nodes)
        names(nodes)[1] <- "Id"
        nodes <- lapply(split(nodes, nodes$Id), function(x){
          x2 <- st_union(x) %>% st_sf(.)
          x2$Id <- unique(x$Id)
          return(x2)}) %>% do.call(rbind, .)
      }
    }
    names(nodes)[1] <- "IdTemp"
    id = "IdTemp"
  }

  if (type %in%  c("centroid", "edge")){
    distance <- euclidean_distances(x = nodes, id = id,
                                    centroid = if(type == "centroid"){TRUE} else {FALSE},
                                    distance_unit = distance_unit,
                                    keep = keep, threshold = threshold,
                                    pairwise = pairwise,
                                    write_table = write)

  } else if (type %in%  c("least-cost", "commute-time")){
    distance <- cost_distances(x = nodes, id = id, LCD = if(type == "least-cost"){TRUE} else {FALSE},
                               resistance = resistance,
                               CostFun = CostFun, ngh = ngh, bounding_circles = bounding_circles ,
                               threshold = threshold, mask = mask, geometry_out = geometry_out,
                               pairwise = pairwise,
                               least_cost.java = least_cost.java,
                               cores.java = cores.java, ram.java = ram.java,
                               write_table = write)

  } else {
    stop("Error, you have to choose a type option")
  }

  return(distance)
}


