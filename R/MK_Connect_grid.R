#' Connectivity indexes in a regular grid
#'
#' Use the function to compute the Protected Connected (ProtConn), EC, PC or IIC indexes in a regular grid.
#' @param nodes \code{sf, SpatVector, SpatialPolygonsDataFrame}. Object containing nodes (e.g., habitat patches or fragments) of each time to analyze information. Nodes are spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.
#' @param area_unit \code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.
#' @param region object of class \code{sf}, \code{SpatialPolygonsDataFrame}. Polygon delimiting the region or study area. It must be
#'  in a projected coordinate system.
#' @param grid \code{list} or object of class \code{sf}, \code{SpatialPolygonsDataFrame}.
#' Use this parameter to generate a grid indicating its characteristics in a \code{list} (see \link[Makurhini]{get_grid}) or enter the name of an sf class \code{sf} or \code{SpatialPolygonsDataFrame} with the grid whose coordinate system must be the same as that of the \code{nodes}.
#'  Example for generating 100 km\out{<sup>2</sup>} hexagons:\cr
#' \code{list(hexagonal = TRUE, cellsize = unit_convert(100, "km2", "m2"), grid_boundary = FALSE, clip = FALSE, tolerance = NULL)}.
#' @param protconn \code{logical}. If \code{TRUE} then the \bold{ProtConn} will be estimated; otherwise, the \bold{PC} index will be estimated.
#' @param distance A \code{list} of parameters to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. \cr
#'  This list must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. For example: \code{distance(type = "least-cost", resistance = raster_resistance)}. \cr
#' To see more arguments see the \link[Makurhini]{distancefile} function.
#' @param distance_threshold A \code{numeric} indicating the dispersal distance (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range.
#' @param probability A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the \code{"PC"} metric. If \code{probability = NULL}, then a probability of 0.5 will be used.
#' @param transboundary \code{numeric}. Buffer to select transboundary polygones, see \link[Makurhini]{MK_ProtConn}.
#' @param intern \code{logical}. Show the progress of the process, \code{default = TRUE}. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @param parallel \code{numeric}. Specify the number of cores to use for parallel processing, \code{default = NULL}. Parallelize the function using furrr package.
#' @references
#' Matt Strimas-Mackey. \url{http://strimas.com/spatial/hexagonal-grids/}.\cr
#' Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the
#' world's ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S. & Torne, J. (2012). Conefor 2.6. Universidad Politecnica de Madrid. Available
#'  at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. (2006). Comparison and development of new graph-based landscape
#'  connectivity indices: towards the priorization of habitat patches and corridors for conservation.
#'  Landscape Ecology, 21(7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity
#' in landscape conservation planning: comparison with existing indices and application to a case study.
#' Landscape and Urban Planning, 83(2-3): 91-103.
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(sf)
#'
#' data("Protected_areas", package = "Makurhini")
#' data("regions", package = "Makurhini")
#' ecoregion <- regions[2,]
#' plot(st_geometry(ecoregion), col="blue")
#' #ProtConn
#' hexagons_priority <- MK_Connect_grid(nodes = Protected_areas,
#'                                     region = ecoregion,
#'                                     area_unit = "ha",
#'                                     grid = list(hexagonal = TRUE,
#'                                                 cellsize = unit_convert(1000, "km2", "m2")),
#'                                     protconn = TRUE,
#'                                     distance_threshold = 3000,
#'                                     probability = 0.5,
#'                                     transboundary = 6000,
#'                                     distance = list(type = "centroid"),
#'                                     intern = TRUE,
#'                                     parallel = NULL)
#' hexagons_priority
#' plot(hexagons_priority["ProtConn"])
#'}
#' @export
#' @importFrom magrittr %>%
#' @importFrom raster crop raster
#' @importFrom sf st_as_sf st_zm st_cast st_buffer st_area st_convex_hull
#' @importFrom future plan multicore multisession availableCores
#' @importFrom furrr future_map_dfr
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom purrr map_df
#' @importFrom methods as

MK_Connect_grid <- function(nodes,
                            area_unit = "ha",
                            region = NULL,
                            grid = list(hexagonal = TRUE, cellsize = NULL,
                                        grid_boundary = FALSE, clip = FALSE, tolerance = NULL),
                            protconn = TRUE,
                            distance_threshold = NULL,
                            probability = NULL,
                            transboundary = NULL,
                            distance = list(type = "centroid"),
                            intern = TRUE, parallel = NULL){
  options(warn = -1)
  . = NULL

  if(!is.null(parallel)){
    if(!is.numeric(parallel)){
      stop("if you use parallel argument then you need a numeric value")
    }
  }

  if(isFALSE(parallel)){
    parallel <- NULL
  }

  if(isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
  }

  message("Step 1. Reviewing parameters")
  base_param1 <- input_grid(node = nodes, landscape = region, unit = area_unit,
                            bdist = if(is.null(transboundary)){0} else{transboundary})

  if(class(base_param1)[1] != "input_grid"){
    stop("error in nodes or region shapefile")
  }

  base_param2 <- metric_class(metric = if(isTRUE(protconn)){"ProtConn"} else {"PC"},
                              distance_threshold = distance_threshold,
                              probability = probability,
                              transboundary = transboundary,
                              distance = distance)

  if(class(base_param2)[1] != "MK_Metric"){
    stop("error in metric parameters")
  } else {
    message("Step 2. Grid processing")
  }

  if(class(grid)[1] == "list"){
    message("Step 2. Grid processing")
    base_param3 <- get_grid(region = region,
                            hexagonal = grid$hexagonal,
                            cellsize = grid$cellsize,
                            grid_boundary = grid$grid_boundary,
                            clip = grid$clip,
                            tolerance = grid$tolerance)
  } else {
    if(any(class(grid)[1] == "sf" | class(grid)[1] == "SpatialPolygonsDataFrame")){
      base_param3 <- get_grid(grid_pol = grid)
    } else {
      stop("You need the grid parameter")
    }
  }

  if(class(base_param3)[1] != "grid"){
    stop("error making the grid")
  }


  base_param4 <- list(base_param1, base_param2, base_param3)

  if(nrow(base_param4[[1]]@nodes) == 0){
    message("Warning message: No nodes found in the region")
  }

  loop <- 1:nrow(base_param4[[3]]@grid)

  if (isTRUE(protconn)) {
    if (isTRUE(intern)) {
      message("Step 3. Processing ProtConn metric on the grid. Progress estimated:")
      pb <- txtProgressBar(0,length(loop), style = 3)
    } else {
      message("Step 3. Processing ProtConn metric on the grid")
    }

    if (!is.null(parallel)){
      works <- as.numeric(availableCores())-1; works <-  if(parallel > works){works}else{parallel}
      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }
      plan(strategy = strat, gc = TRUE, workers = works)
      result_1 <- tryCatch(future_map_dfr(loop, function(x){

        LA <- st_area(base_param4[[3]]@grid[x,]) %>%
          unit_convert(., "m2", base_param4[[1]]@area_unit)

        #nodes and distances,
        nodes.1 <- tryCatch(Protconn_nodes(x = base_param4[[3]]@grid[x,],
                                           y = base_param4[[1]]@nodes,
                                           buff = base_param4[[2]]@transboundary,
                                           xsimplify = FALSE,
                                           metrunit = base_param4[[1]]@area_unit,
                                           protconn_bound = FALSE), error = function(err)err)

        if(inherits(nodes.1, "error")){
          stop(paste0("error first nodes selection. Check grid: ", x))
        }
        #Tiene 2 o más nodos dentro de la region
        if(is.list(nodes.1)){
          if(base_param4[[2]]@distance$type %in% c("least-cost", "commute-time")){
            if(is.null(base_param4[[2]]@distance$resistance)){
              stop("error, you need a resistance raster")
            } else {
              centroid <- st_centroid(nodes.1[[1]])
              mask <- st_convex_hull(st_union(centroid)) %>%
                st_buffer(res(base_param4[[2]]@distance$resistance)[1]*30)
              resist <- crop(base_param4[[2]]@distance$resistance, as(mask, 'Spatial'))
            }
          } else {
            resist <- NULL
          }

          distance.1 <- tryCatch(protconn_dist(nodes.1[[1]], id = "OBJECTID",
                                               y = base_param4[[2]]@distance,
                                               r = base_param4[[1]]@region,
                                               resistance = resist),
                                 error = function(err)err)
          if(inherits(distance.1, "error")){
            stop(paste0("error distance. Check", " grid: ", x))
          }

          ProtConn_grid <- get_protconn_grid(x = nodes.1,
                                             y = distance.1,
                                             p = base_param4[[2]]@probability,
                                             pmedian = TRUE,
                                             d = base_param4[[2]]@distance_threshold,
                                             LA = LA, bound = FALSE)
          ProtConn_grid <- round(ProtConn_grid, 5)

        } else if(is.numeric(nodes.1)){
          ProtConn_grid <- data.frame(ECA = if(nodes.1 > LA){LA}else{nodes.1},
                                      PC = if(nodes.1 >= LA){1}else{nodes.1/LA^2},
                                      LA = LA,
                                      Protected.surface = nodes.1,
                                      Prot = if((100 * (nodes.1 / LA)) > 100){100}else{100 * (nodes.1/LA)},
                                      Unprotected = if((100 - (100 * (nodes.1 / LA))) < 0){0}else{100 - (100 * (nodes.1 / LA))},
                                      ProtConn = if((100 * (nodes.1 / LA)) > 100){100}else{100 * (nodes.1 / LA)},
                                      ProtUnconn = 0,
                                      RelConn = NA,
                                      ProtConn_Prot = 100,
                                      ProtConn_Trans = NA,
                                      ProtConn_Unprot = NA,
                                      ProtConn_Within = 100,
                                      ProtConn_Contig = NA,
                                      ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                      ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)
          ProtConn_grid[,which(!is.na(ProtConn_grid))] <- round(ProtConn_grid[,which(!is.na(ProtConn_grid))], 5)
        } else {
          ProtConn_grid <- data.frame(ECA = NA,
                                      PC = NA,
                                      LA = LA,
                                      Protected.surface = 0,
                                      Prot = 0,
                                      Unprotected = 100,
                                      ProtConn = NA,
                                      ProtUnconn = NA,
                                      RelConn = NA,
                                      ProtConn_Prot = NA,
                                      ProtConn_Trans = NA,
                                      ProtConn_Unprot = NA,
                                      ProtConn_Within = NA,
                                      ProtConn_Contig = NA,
                                      ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                      ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)
        }

        return(ProtConn_grid) }, .progress = intern), error = function(err) err)
      close_multiprocess(works)
    } else {
      result_1 <- tryCatch(map_df(loop, function(x) {
        if (isTRUE(intern)) {
          setTxtProgressBar(pb, x)
        }
        LA <- st_area(base_param4[[3]]@grid[x,]) %>%
          unit_convert(., "m2", base_param4[[1]]@area_unit)

        #nodes and distances,
        nodes.1 <- tryCatch(Protconn_nodes(x = base_param4[[3]]@grid[x,],
                                           y = base_param4[[1]]@nodes,
                                           buff = base_param4[[2]]@transboundary,
                                           xsimplify = FALSE,
                                           metrunit = base_param4[[1]]@area_unit,
                                           protconn = TRUE,
                                           protconn_bound = FALSE), error = function(err)err)

        if(inherits(nodes.1, "error")){
          stop(paste0("error first nodes selection. Check grid: ", x))
        }

        if(is.list(nodes.1)){
          if(base_param4[[2]]@distance$type %in% c("least-cost", "commute-time")){
            if(is.null(base_param4[[2]]@distance$resistance)){
              stop("error, you need a resistance raster")
            } else {
              centroid <- st_centroid(nodes.1[[1]])
              mask <- st_convex_hull(st_union(centroid)) %>%
                st_buffer(res(base_param4[[2]]@distance$resistance)[1]*30)
              resist <- crop(base_param4[[2]]@distance$resistance, as(mask, 'Spatial'))
            }
          } else {
            resist <- NULL
          }

          distance.1 <- tryCatch(protconn_dist(nodes.1[[1]], id = "OBJECTID",
                                               y = base_param4[[2]]@distance,
                                               r = base_param4[[1]]@region,
                                               resistance = resist),
                                 error = function(err)err)
          if(inherits(distance.1, "error")){
            stop(paste0("error distance. Check", " grid: ", x))
          }

          ProtConn_grid <- get_protconn_grid(x = nodes.1,
                                             y = distance.1,
                                             p = base_param4[[2]]@probability,
                                             pmedian = TRUE,
                                             d = base_param4[[2]]@distance_threshold,
                                             LA = LA, bound = FALSE)
          ProtConn_grid <- round(ProtConn_grid, 5)

        } else if(is.numeric(nodes.1)){
          ProtConn_grid <- data.frame(ECA = if(nodes.1 >= LA){LA}else{nodes.1},
                                      PC = if(nodes.1 >= LA){1}else{nodes.1/LA^2},
                                      LA = LA,
                                      Protected.surface = nodes.1,
                                      Prot = if((100 * (nodes.1 / LA)) > 100){100}else{100 * (nodes.1/LA)},
                                      Unprotected = if((100 - (100 * (nodes.1 / LA))) < 0){0}else{100 - (100 * (nodes.1 / LA))},
                                      ProtConn = if((100 * (nodes.1 / LA)) > 100){100}else{100 * (nodes.1 / LA)},
                                      ProtUnconn = 0,
                                      RelConn = NA,
                                      ProtConn_Prot = 100,
                                      ProtConn_Trans = NA,
                                      ProtConn_Unprot = NA,
                                      ProtConn_Within = 100,
                                      ProtConn_Contig = NA,
                                      ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                      ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)
          ProtConn_grid[,which(!is.na(ProtConn_grid))] <- round(ProtConn_grid[,which(!is.na(ProtConn_grid))], 5)

        } else {
          ProtConn_grid <- data.frame(ECA = NA,
                                      PC = NA,
                                      LA = LA,
                                      Protected.surface = 0,
                                      Prot = 0,
                                      Unprotected = 100,
                                      ProtConn = NA,
                                      ProtUnconn = NA,
                                      RelConn = NA,
                                      ProtConn_Prot = NA,
                                      ProtConn_Trans = NA,
                                      ProtConn_Unprot = NA,
                                      ProtConn_Within = NA,
                                      ProtConn_Contig = NA,
                                      ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                      ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)
        }
        return(ProtConn_grid)
      }), error = function(err) err)
    }
  } else {
    if (isTRUE(intern)) {
      message("Step 3. Processing PC and ECA metrics on the grid. Progress estimated:")
      pb <- txtProgressBar(0,length(loop), style = 3)
    } else {
      message("Step 3. Processing PC and ECA metrics on the grid")
    }

    if (!is.null(parallel)) {
      works <- as.numeric(availableCores())-1; works <-  if(parallel > works){works}else{parallel}
      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }
      plan(strategy = strat, gc = TRUE, workers = works)
      result_1 <- tryCatch(future_map_dfr(loop, function(x) {
        LA <- st_area(base_param4[[3]]@grid[x,])
        LA <- unit_convert(LA, "m2", base_param4[[1]]@area_unit)

        #nodes and distances,
        nodes.1 <- tryCatch(Protconn_nodes(x = base_param4[[3]]@grid[x,],
                                           y = base_param4[[1]]@nodes,
                                           buff = base_param4[[2]]@transboundary,
                                           xsimplify = FALSE,
                                           metrunit = base_param4[[1]]@area_unit,
                                           protconn = FALSE,
                                           protconn_bound = FALSE), error = function(err)err)

        if(inherits(nodes.1, "error")){
          stop(paste0("error first nodes selection. Check grid: ", x))
        }

        if(is.list(nodes.1)){
          if(base_param4[[2]]@distance$type %in% c("least-cost", "commute-time")){
            if(is.null(base_param4[[2]]@distance$resistance)){
              stop("error, you need a resistance raster")
            } else {
              centroid <- st_centroid(nodes.1[[1]])
              mask <- st_convex_hull(st_union(centroid)) %>%
                st_buffer(res(base_param4[[2]]@distance$resistance)[1]*30)
              resist <- crop(base_param4[[2]]@distance$resistance, as(mask, 'Spatial'))
            }
          } else {
            resist <- NULL
          }

          distance.1 <- tryCatch(protconn_dist(nodes.1[[1]], id = "OBJECTID",
                                               y = base_param4[[2]]@distance,
                                               r = base_param4[[1]]@region,
                                               resistance = resist),
                                 error = function(err)err)
          if(inherits(distance.1, "error")){
            stop(paste0("error distance. Check", " grid: ", x))
          }

          PC_grid <- get_pc_grid(x = nodes.1[[1]],
                                 y = distance.1,
                                 p = base_param4[[2]]@probability,
                                 pmedian = TRUE,
                                 d = base_param4[[2]]@distance_threshold,
                                 LA = LA) %>% round(., 5)

        } else if(is.numeric(nodes.1)){
          PC_grid <- data.frame(Protected.surface = nodes.1,
                                LA = LA,
                                ECA = if(nodes.1 >= LA){LA}else{nodes.1},
                                ECA.Normalized = if(nodes.1 >= LA){100}else{(nodes.1*100)/LA},
                                PC = if(nodes.1 >= LA){1}else{nodes.1/LA^2})
          PC_grid[,c(1:4)] <- round(PC_grid[,c(1:4)], 5)
        } else {
          PC_grid <- data.frame(Protected.surface = 0,
                                LA = LA,
                                ECA = NA,
                                ECA.Normalized = NA,
                                PC = NA)
          PC_grid[,c(1:2)] <- round(PC_grid[,c(1:2)], 5)
        }

        return(PC_grid)}, .progress = intern), error = function(err) err)
      close_multiprocess(works)
    } else {
      result_1 <- tryCatch(map_df(loop, function(x) {
        if (isTRUE(intern)) {
          setTxtProgressBar(pb, x)
        }

        LA <- st_area(base_param4[[3]]@grid[x,]) %>%
          unit_convert(., "m2", base_param4[[1]]@area_unit)

        #nodes and distances,
        nodes.1 <- tryCatch(Protconn_nodes(x = base_param4[[3]]@grid[x,],
                                           y = base_param4[[1]]@nodes,
                                           buff = NULL,
                                           xsimplify = FALSE,
                                           metrunit = base_param4[[1]]@area_unit,
                                           protconn = FALSE,
                                           protconn_bound = FALSE), error = function(err)err)

        if(inherits(nodes.1, "error")){
          stop(paste0("error first nodes selection. Check grid: ", x))
        }

        if(is.list(nodes.1)){
          if(base_param4[[2]]@distance$type %in% c("least-cost", "commute-time")){
            if(is.null(base_param4[[2]]@distance$resistance)){
              stop("error, you need a resistance raster")
            } else {
              centroid <- st_centroid(nodes.1[[1]])
              mask <- st_convex_hull(st_union(centroid)) %>%
                st_buffer(res(base_param4[[2]]@distance$resistance)[1]*30)
              resist <- crop(base_param4[[2]]@distance$resistance, as(mask, 'Spatial'))
            }
          } else {
            resist <- NULL
          }

          distance.1 <- tryCatch(protconn_dist(nodes.1[[1]], id = "OBJECTID",
                                               y = base_param4[[2]]@distance,
                                               r = base_param4[[1]]@region,
                                               resistance = resist),
                                 error = function(err)err)
          if(inherits(distance.1, "error")){
            stop(paste0("error distance. Check", " grid: ", x))
          }

          PC_grid <- get_pc_grid(x = nodes.1[[1]],
                                 y = distance.1,
                                 p = base_param4[[2]]@probability,
                                 pmedian = TRUE,
                                 d = base_param4[[2]]@distance_threshold,
                                 LA = LA) %>% round(., 5)

        } else if(is.numeric(nodes.1)){
          PC_grid <- data.frame(Protected.surface = nodes.1,
                                LA = LA,
                                ECA = if(nodes.1 >= LA){LA}else{nodes.1},
                                ECA.Normalized = if(nodes.1 >= LA){100}else{(nodes.1*100)/LA},
                                PC = if(nodes.1 >= LA){1}else{nodes.1/(LA^2)})
          PC_grid[,c(1:4)] <- round(PC_grid[,c(1:4)], 5)
        } else {
          PC_grid <- data.frame(Protected.surface = 0,
                                LA = LA,
                                ECA = NA,
                                ECA.Normalized = NA,
                                PC = NA)
          PC_grid[,c(1:2)] <- round(PC_grid[,c(1:2)], 5)
        }

        return(PC_grid)}), error = function(err) err)
    }
  }

  if(inherits(result_1, "error")){
    stop("error, check your input files, there may be no nodes in the region")
  }

  result_2 <- cbind(base_param4[[3]]@grid, result_1); result_2$IdTemp  <- NULL
  return(result_2)
}
