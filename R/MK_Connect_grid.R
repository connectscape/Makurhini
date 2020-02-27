#' Connectivity indexes in a regular grid
#'
#' Use the function to compute the Protected Connected (ProtConn)or EC, PC, IIC indexes on a grid.
#' @param nodes	object of class sf, sfc, sfg or SpatialPolygons. Nodes shapefile, the shapefile must be
#'  in a projected coordinate system.
#' @param region object of class sf, sfc, sfg or SpatialPolygons. Region shapefile, the shapefile must be
#'  in a projected coordinate system.
#' @param grid_pol object of class sf, sfc, sfg or SpatialPolygons. Grid hexagones or squares. The
#' shapefile must be in a projected coordinate system.
#' @param grid_id character. Column name of the grid ID.
#' @param grid_type character. If grid is null you can make a regular grid of "hexagonal" or "square".
#' @param cellsize numeric. Grid area (square kilometers).
#' @param grid_boundary logical.If TRUE, the Incomplete "hexagons" or "squares" in the boundaries of
#' the region will be discarded
#' @param clip logical. If TRUE, the new grid will be clipped to the region area. The operation time
#' will be longer the greater the number of vertices in the polygon of the region,
#' if it is a region with many vertices use the argument "tolerance".
#' @param tolerance numeric. If "clip" is equal to TRUE reduces the number of vertices in the region polygon.
#' @param distance list. See distancefile(). Example, list(type= "centroid", resistance = NULL).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5
#' (default) that is 50 percentage of probability connection. Use in case of selecting the "PC"
#' metric or "ProtConn".
#' @param distance_thresholds numeric. Distance thresholds to establish connections (meters).
#' @param ProtConn list. ProtConn parameters (see, "Makurhini::MK_ProtConn"). Example,
#' list(attribute = "Intersected area", transboundary = NULL).
#' @param PC_IIC list. PC or IIC metric parameters (see, "Makurhini::MK_dPCIIC"). Example, list(attribute  = NULL,
#' area_unit = "ha", metric = c("IIC", "PC")).
#' @param intern logical. Show the progress of the process, default = TRUE.
#' @param parallel logical. Parallelize the function using furrr package and multiprocess
#' plan, default = FALSE.
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
#' @export
#' @examples
#' \dontrun{
#' data("Protected_areas", package = "Makurhini")
#' #plot(Protected_areas, col="green")
#'
#' data("regions", package = "Makurhini")
#' region <- regions[2,]
#' plot(region, col="blue")
#'
#' hexagons_priority <- MK_Connect_grid(nodes = Protected_areas, region = region,
#'                                   grid_type = "hexagonal", cellsize = 3000,
#'                                   distance = list(type = "centroid"),
#'                                   distance_thresholds = 10000, clip = FALSE,
#'                                   PC_IIC = list(metric = "PC"),
#'                                   intern = TRUE, parallel = TRUE)
#' hexagons_priority
#'
#'
#' hexagons_priority_2 <- MK_Connect_grid(nodes = Protected_areas, region = region,
#'                                    grid_type = "hexagonal", cellsize = 3000,
#'                                    distance = list(type = "centroid"),
#'                                    distance_thresholds = 5000, clip = FALSE,
#'                                    ProtConn = list(attribute ="Intersected area",
#'                                                    transboundary = 10000),
#'                                    intern = TRUE, parallel = TRUE)
#' hexagons_priority_2
#' }
#' @importFrom magrittr %>%
#' @importFrom raster extent projection projection<- raster
#' @import sf
#' @importFrom rmapshaper ms_dissolve
#' @importFrom sp spsample HexPoints2SpatialPolygons
#' @importFrom future plan multiprocess
#' @importFrom furrr future_map
#' @importFrom dplyr progress_estimated
#' @importFrom purrr map

MK_Connect_grid <- function(nodes, region = NULL, grid_pol = NULL, grid_id = NULL,
                         grid_type = c("hexagonal", "square"),
                         cellsize = NULL, grid_boundary = FALSE,
                         clip = FALSE, tolerance = NULL,
                         distance = list(type = "centroid"),
                         probability = 0.5,
                         distance_thresholds = NULL,
                         ProtConn = list(attribute = "Intersected area",
                                         area_unit = "ha",
                                         transboundary = NULL),
                         PC_IIC = list(attribute  = NULL,
                                       area_unit = "ha",
                                       metric = "PC"),
                         intern = TRUE, parallel = FALSE){
  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (is.null(region)) {
    if (is.null(grid_pol)) {
      stop("error missing shapefile file of region or grid_pol")
    }
  } else {
    if (is.numeric(region) | is.character(region)) {
      stop("error missing shapefile file of region")
    }

    if (is.null(grid_pol)) {
      if (is.null(cellsize)) {
        stop("error missing cellsize(km2)")
      }
      if (!grid_type %in% c("hexagonal", "square")) {
        stop("Type must be either 'hexagonal' or 'square'")
      }
    }
  }

  if (is.null(distance_thresholds)) {
    stop("error missing numeric distance threshold(s)")
  }

  if (missing(PC_IIC)) {
    area_unit <- ProtConn$area_unit
    if (is.null(probability) | !is.numeric(probability)) {
      stop("error missing probability")
    }
  } else {
    area_unit <- PC_IIC$area_unit
    if (!PC_IIC$metric %in% c("PC", "IIC")) {
      stop("Type must be either 'PC', or 'IIC'")
    }
    if (PC_IIC$metric == "PC") {
      if (is.null(probability) | !is.numeric(probability)) {
        stop("error missing probability")
      }
    }
  }

  if(is.null(area_unit)){
    area_unit = "ha"
  }


  options(warn = -1)
  ttt.2 <- getwd()

  #number of vertices
  num_vert <- function(y){
    if(class(y)[1] == "SpatialPolygonsDataFrame") {
      y <- st_as_sf(y)
    }
    y <- st_cast(y$geometry, "MULTIPOINT")
    y <- sapply(y, length)
    y <- sum(y)/2
    return(y)
  }

  make_grid <- function(x, type, cell_width, cell_area,
                        clip = FALSE, tolerance = tolerance, grid_boundary = grid_boundary) {
    if (class(x)[1] == "sf") {
      x <- st_zm(x)
      x <- as(x, "Spatial")
    }

    if (!type %in% c("square", "hexagonal")) {
      stop("Type must be either 'square' or 'hexagonal'")
    }

    if (missing(cell_width)) {
      if (missing(cell_area)) {
        stop("Must provide cell_width or cell_area")
      } else {
        if (type == "square") {
          cell_width <- sqrt(cell_area)
        } else if (type == "hexagonal") {
          cell_width <- sqrt(2 * cell_area/sqrt(3))
        }
      }
    }
    ext <- as(extent(x) + cell_width, "SpatialPolygons")
    projection(ext) <- projection(x)

    if (type == "square") {
      g <- raster(ext, resolution = cell_width)
      g <- as(g, "SpatialPolygons") %>% st_as_sf() %>%
         st_cast("POLYGON")
    } else if (type == "hexagonal") {
      g <- spsample(ext, type = "hexagonal", cellsize = cell_width,
                    offset = c(0, 0))
      g <- HexPoints2SpatialPolygons(g, dx = cell_width) %>%
        st_as_sf() %>% st_cast("POLYGON")
    } else {
      stop("Select a grid_type = c('square', 'hexagonal')")
    }

    x2 <- st_as_sf(x) %>% st_cast("POLYGON")

    if (num_vert(x2) > 100) {
      region_1 <- tryCatch(st_simplify(x2, dTolerance = ((extent(x)[4]-extent(x)[3])*2)/100,
                                       preserveTopology = TRUE), error = function(err) err)

      if (inherits(region_1, "error")) {
        x2 <- st_buffer(x2, dist = 0)
        x2 <- st_simplify(x2, dTolerance = ((extent(x)[4]-extent(x)[3])*2)/100, preserveTopology = TRUE)
        } else {
        x2 <- region_1
      }
    }


    x2 <- tryCatch(x2 %>% st_cast("POLYGON") %>% ms_dissolve() %>%
                     st_buffer(dist = 0), error = function(err) err)

    if (inherits(x2, "error")) {
      x2 <- tryCatch(x2 %>% st_buffer(., dist = 0) %>%
                       ms_dissolve() %>% st_buffer(., dist = 0), error = function(err) err)
    }

    if (isTRUE(clip)) {
      if(is.null(tolerance)){
        if (num_vert(x2) > 100) {
          region_1 <- tryCatch(st_simplify(x2, dTolerance = ((extent(x)[4]-extent(x)[3])*0.1)/100,
                                           preserveTopology = TRUE), error = function(err) err)

          if (inherits(region_1, "error")) {
            x2 <- st_buffer(x2, dist = 0)
            x2 <- st_simplify(x2, dTolerance =((extent(x)[4]-extent(x)[3])*0.1)/10, preserveTopology = TRUE)
          } else {
            x2 <- region_1
          }

        }
      } else {
        region_1 <- tryCatch(st_simplify(x2, dTolerance = tolerance,
                                         preserveTopology = TRUE), error = function(err) err)
        if (inherits(region_1, "error")) {
          x2 <- st_buffer(x2, dist = 0)
          x2 <- st_simplify(x2, dTolerance = tolerance, preserveTopology = TRUE)
        } else {
          x2 <- region_1
        }
      }

      g2 <- tryCatch(st_intersection(g, y = x2) %>% st_cast("POLYGON"),
                     error = function(err) err)

      if (inherits(g2, "error")) {
        g2 <- tryCatch(st_intersection(g, y = st_buffer(x2,dist = 0)) %>%
                         st_cast("POLYGON"), error = function(err) err)
      }

    } else {
      a <- unique(format(round(as.numeric(st_area(g)), 0), scientific = F))

      g2 <- MK_selectbyloc(g, x2)

      if(isTRUE(grid_boundary)){
        if (num_vert(x2) > 100) {
          region_1 <- tryCatch(st_simplify(x2, dTolerance = ((extent(x)[4]-extent(x)[3])*2)/100,
                                           preserveTopology = TRUE), error = function(err) err)
          if (inherits(region_1, "error")) {
            x2 <- st_buffer(x2, dist = 0)
            x2 <- st_simplify(x2, dTolerance = ((extent(x)[4]-extent(x)[3])*2)/100, preserveTopology = TRUE)
          } else {
            x2 <- region_1
          }
        }

        g3 <- tryCatch(st_intersection(g2, y = x2) %>% st_cast("POLYGON"),
                       error = function(err) err)

        if (inherits(g3, "error")) {
          g3 <- tryCatch(st_intersection(g2, y = st_buffer(x2,dist = 0)) %>%
                          st_cast("POLYGON"), error = function(err) err)
        }

        g3$area <- format(round(as.numeric(st_area(g3)),0), scientific = F)
        g3 <- g3[which(g3$area == a), ]
        g2 <- g3
      }
    }

    row.names(g2) <- as.character(1:nrow(g2))
    g2$rmapshaperid <- NULL
    g2$area <- NULL
    return(g2)
  }

  ###
  if (class(nodes)[1] == "SpatialPolygonsDataFrame") {
    nodes <- st_as_sf(x = nodes) %>% st_zm()
  } else {
    nodes <- st_zm(nodes)
  }



  if (is.null(grid_pol) & !is.null(region)) {
    cellsize <- cellsize*1000000
    x_grid <- make_grid(x = region, type = grid_type, cell_area = cellsize, clip = clip, tolerance = tolerance, grid_boundary = grid_boundary)
  } else {
    if (class(grid_pol)[1] == "SpatialPolygonsDataFrame") {
      x_grid <- st_as_sf(grid_pol) %>% st_zm() %>% st_cast("POLYGON")
    } else {
      x_grid <- grid_pol %>% st_zm() %>% st_cast("POLYGON")
    }
  }

  ###grid to list
  x_grid$TempID <- 1:nrow(x_grid)
  x_grid <- x_grid["TempID"]
  x_grid <- split(x_grid, x_grid$TempID)
  if (missing(PC_IIC) & !missing(ProtConn)) {
    if (isTRUE(parallel)) {
      plan(strategy = multiprocess)
      resultado_1 <- tryCatch(future_map(x_grid, function(x) {
        test1 <- MK_ProtConn(nodes = nodes, region = x,
                             area_unit = ProtConn$area_unit,
                             res_attribute = ProtConn$res_attribute,
                             thintersect = ProtConn$thintersect, attribute = ProtConn$attribute,
                             distance = distance, distance_thresholds = distance_thresholds,
                             probability = probability, transboundary = ProtConn$transboundary,
                             LA = NULL, plot = FALSE, dPC = FALSE, write = NULL,
                             intern = FALSE)
        test2 <- test1[[1]]
        n <- as.vector(test2[[1]][[3]])
        x2 <- as.data.frame(test2[[1]][[4]])
        x3 <- as.data.frame(t(x2))
        names(x3) <- n
        x3$TempID <- paste0(x[["TempID"]])
        return(x3)
      }, .progress = intern), error = function(err) err)
    } else {
      pb <- progress_estimated(length(x_grid), 0)
      resultado_1 <- tryCatch(map(x_grid, function(x) {
        if (isTRUE(intern)) {
          pb$tick()$print()
        }
        test1 <- MK_ProtConn(nodes = nodes, region = x,
                             area_unit = ProtConn$area_unit,
                             res_attribute = ProtConn$res_attribute,
                             thintersect = ProtConn$thintersect, attribute = ProtConn$attribute,
                             distance = distance, distance_thresholds = distance_thresholds,
                             probability = probability, transboundary = ProtConn$transboundary,
                             LA = NULL, plot = FALSE, dPC = FALSE, write = NULL,
                             intern = FALSE)
        test2 <- test1[[1]]
        n <- as.vector(test2[[1]][[3]])
        x2 <- as.data.frame(test2[[1]][[4]])
        x3 <- as.data.frame(t(x2))
        names(x3) <- n
        x3$TempID <- paste0(x[["TempID"]])
        return(x3)
      }), error = function(err) err)
    }
  } else if (!missing(PC_IIC) & missing(ProtConn)) {
    if (isTRUE(parallel)) {
     nodes1 <- nodes
     nodes1 <- st_cast(nodes1, "POLYGON")
     plan(strategy = multiprocess)
      resultado_1 <- tryCatch(future_map(x_grid, function(x) {
        area_x <- unit_convert(as.numeric(st_area(x)), "m2", area_unit)
        if (is.null(PC_IIC$attribute)) {
          nodes2 <- MK_selectbyloc(nodes1, x, id = NULL, area_unit = area_unit,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]
          if (nrow(nodes2) > 1) {
            nodes3 <- st_intersection(nodes2, x) %>%
              ms_dissolve(.) %>% st_buffer(., dist = 0)
            nodes2_r <- nodes2
            nodes2 <- st_cast(nodes3, "POLYGON") %>%
              st_buffer(., dist = 0)
            nodes2$attrib <- unit_convert(as.numeric(st_area(nodes2)), "m2", area_unit)
            if (nrow(nodes2) == 1) {
              nodes2 <- (nodes2$attrib * 100)/area_x
            } else if (nrow(nodes2) == 0) {
              nodes2 = NULL
            }

          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- (nodes2$Area2 * 100)/area_x
          } else {
            nodes2 = NULL
          }
        } else {
          nodes2 <- MK_selectbyloc(nodes1, x, id = NULL, area_unit = area_unit,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]

          if (nrow(nodes2) > 1) {
            nodes3 <- st_intersection(nodes2, x) %>%
              ms_dissolve(., field = PC_IIC$attribute) %>%
              st_buffer(., dist = 0)

            nodes3 <- st_cast(nodes3, "POLYGON") %>%
              st_buffer(., dist = 0)

            if (nrow(nodes3) > 1) {
              nodes3$IDTemp <- 1:nrow(nodes3)
              d <- distancefile(nodes3, id = "IDTemp", type = "edge")

              if (unique(d$Distance) == 0) {
                nodes2_r <- nodes2
                nodes2 <- unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes2 <- nodes3[ PC_IIC$attribute][[1]] * nodes2
                }
                nodes2 <- sum((nodes2 * 100)/area_x)
              } else {
                nodes3$attrib <- nodes3[PC_IIC$attribute][[1]]
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes3$attrib <- nodes3$attrib * unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
                }
                nodes2 <- nodes3
              }
            } else {
              nodes2_r <- nodes2
              nodes2 <-unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
              if (isTRUE(PC_IIC$attribute_weighted_area)) {
                nodes2 <- nodes3[ PC_IIC$attribute][[1]] * nodes2
              }
              nodes2 <- sum((nodes2 * 100)/area_x)
            }
          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- nodes2$PercPi
            nodes3 <- nodes2[PC_IIC$attribute][[1]]
            if (isTRUE(PC_IIC$attribute_weighted_area)) {
              nodes3 <- nodes3 * nodes2$Area2
            }
            nodes2 <- (nodes3 * 100)/area_x
          } else {
            nodes2 = NULL
          }
        }
        ####
        if (class(nodes2)[1] == "numeric" | is.null(nodes2)) {
          if (is.null(nodes2)) {
            x3 <- c(NA, NA, NA) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(PC_IIC$metric, "EC", "Normalized_EC")
            x3$PArea <- 0
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          } else {

            x3 <- c(NA, sum(unit_convert(as.numeric(st_area(nodes2_r)), "m2", area_unit)), 100) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(PC_IIC$metric, "EC", "Normalized_EC")

            x3$PArea <- if(((sum(nodes2_r$Area2 * 100))/area_x) > 100){
              100
            } else {
              sum((nodes2_r$Area2 * 100))/area_x
            }

            x3$EC <- if (sum(nodes2_r$Area2) > area_x) {
              area_x
            } else {
              sum(nodes2_r$Area2)
            }
            x3$Normalized_EC <- x3$PArea
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          }
        } else {
          test <- MK_dPCIIC(nodes = nodes2, id = NULL,
                            attribute = "attrib", restauration = NULL,
                            distance = distance, metric = PC_IIC$metric,
                            probability = probability, distance_thresholds = distance_thresholds,
                            overall = TRUE, LA = area_x,
                            write = NULL)
          test <- test[[1]][[2]]
          x2 <- as.data.frame(test[, 2])
          x3 <- as.data.frame(t(x2))
          x3 <- x3[, c(2:3)]
          names(x3) <- c("EC", PC_IIC$metric)
          x3$Normalized_EC <- (x3$EC * 100)/area_x
          x3$PArea <- (sum(nodes2$attrib) * 100)/area_x
          x3$TempID <- paste0(x[["TempID"]])
          x3 <- x3[, c(4, 2, 1, 3, 5)]
        }
        return(x3)
      }, .progress = intern), error = function(err) err)
    } else {
      pb <- progress_estimated(length(x_grid), 0)
      nodes1 <- nodes
      nodes1 <- st_cast(nodes1, "POLYGON")
      resultado_1 <- tryCatch(map(x_grid, function(x) {
        if (isTRUE(intern)) {
          pb$tick()$print()
        }
        area_x <- unit_convert(as.numeric(st_area(x)), "m2", area_unit)
        if (is.null(PC_IIC$attribute)) {
          nodes2 <- MK_selectbyloc(nodes1, x, id = NULL, area_unit = area_unit,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]
          if (nrow(nodes2) > 1) {
            nodes3 <- st_intersection(nodes2, x) %>%
              ms_dissolve(.) %>% st_buffer(., dist = 0)
            nodes2_r <- nodes2
            nodes2 <- st_cast(nodes3, "POLYGON") %>%
              st_buffer(., dist = 0)
            nodes2$attrib <- unit_convert(as.numeric(st_area(nodes2)), "m2", area_unit)

            if (nrow(nodes2) == 1) {
              nodes2 <- (nodes2$attrib * 100)/area_x
            } else if (nrow(nodes2) == 0) {
              nodes2 = NULL
            }

          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- (nodes2$Area2 * 100)/area_x
          } else {
            nodes2 = NULL
          }
        } else {
          nodes2 <- MK_selectbyloc(nodes1, x, id = NULL, area_unit = area_unit,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]

          if (nrow(nodes2) > 1) {
            nodes3 <- st_intersection(nodes2, x) %>%
              ms_dissolve(., field = PC_IIC$attribute) %>%
              st_buffer(., dist = 0)

            nodes3 <- st_cast(nodes3, "POLYGON") %>%
              st_buffer(., dist = 0)

            if (nrow(nodes3) > 1) {
              nodes3$IDTemp <- 1:nrow(nodes3)
              d <- distancefile(nodes3, id = "IDTemp", type = "edge")

              if (unique(d$Distance) == 0) {
                nodes2_r <- nodes2
                nodes2 <-unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes2 <- nodes3[ PC_IIC$attribute][[1]] * nodes2
                }
                nodes2 <- sum((nodes2 * 100)/area_x)
              } else {
                nodes3$attrib <- nodes3[PC_IIC$attribute][[1]]
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes3$attrib <- nodes3$attrib * unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
                }
                nodes2 <- nodes3
              }
            } else {
              nodes2_r <- nodes2
              nodes2 <- unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
              if (isTRUE(PC_IIC$attribute_weighted_area)) {
                nodes2 <- nodes3[PC_IIC$attribute][[1]] * nodes2
              }
              nodes2 <- sum((nodes2 * 100)/area_x)
            }
          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- nodes2$PercPi
            nodes3 <- nodes2[PC_IIC$attribute][[1]]
            if (isTRUE(PC_IIC$attribute_weighted_area)) {
              nodes3 <- nodes3 * nodes2$Area2
            }
            nodes2 <- (nodes3 * 100)/area_x
          } else {
            nodes2 = NULL
          }
        }
        ##
        if (class(nodes2)[1] == "numeric" | is.null(nodes2)) {
          if (is.null(nodes2)) {
            x3 <- c(NA, NA, NA) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(PC_IIC$metric, "EC", "Normalized_EC")
            x3$PArea <- 0
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          } else {

            x3 <- c(NA, sum(unit_convert(as.numeric(st_area(nodes2_r)), "m2", area_unit)), 100) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(PC_IIC$metric, "EC", "Normalized_EC")
            x3$PArea <- if(((sum(nodes2_r$Area2 * 100))/area_x) > 100){
              100
            } else {
              sum((nodes2_r$Area2 * 100))/area_x
            }

            x3$EC <- if (sum(nodes2_r$Area2) > area_x) {
              area_x
            } else {
              sum(nodes2_r$Area2)
            }
            x3$Normalized_EC <- x3$PArea
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          }
        } else {
          test <- MK_dPCIIC(nodes = nodes2, id = NULL,
                                attribute = "attrib", restauration = NULL,
                                distance = distance, metric = PC_IIC$metric,
                                probability = probability, distance_thresholds = distance_thresholds,
                                overall = TRUE, LA = area_x, write = NULL)
          test <- test[[1]][[2]]
          x2 <- as.data.frame(test[, 2])
          x3 <- as.data.frame(t(x2))
          x3 <- x3[, c(2:3)]
          names(x3) <- c("EC", PC_IIC$metric)
          x3$Normalized_EC <- (x3$EC * 100)/area_x
          x3$PArea <- (sum(nodes2$attrib) * 100)/area_x
          x3$TempID <- paste0(x[["TempID"]])
          x3 <- x3[, c(4, 2, 1, 3, 5)]
        }
        return(x3)
      }), error = function(err) err)
    }
  } else {
    stop("check arguments or select only one argument ProtConn or PC_IIC")
  }
  ####
  if (inherits(resultado_1, "error")) {
    setwd(ttt.2)
  } else {
    result_2 <- do.call(rbind, resultado_1)
    result_2 <- result_2[, c(ncol(result_2), 1:(ncol(result_2) - 1))]
    x_grid2 <- do.call(rbind, x_grid)
    x_grid2 <- merge(x_grid2, result_2, by = "TempID")
    if (!is.null(grid_id)) {
      x_grid2$TempID <- grid_pol[, which(names(grid_pol) == grid_id)][[1]]
      names(x_grid2)[1] <- grid_id
    } else {
      names(x_grid2)[1] <- "id"
    }
    resultado_1 <- x_grid2
    setwd(ttt.2)
  }
  return(resultado_1)
}
