#' Connectivity indexes in a regular grid
#'
#' Use the function to compute the Protected Connected (ProtConn)or EC, PC, IIC indexes on a grid.
#' @param nodes	object of class sf, sfc, sfg or SpatialPolygons. nodes shapefile. The shapefile must be in a projected coordinate system.
#' @param region object of class sf, sfc, sfg or SpatialPolygons. region shapefile.The shapefile must be in a projected coordinate system.
#' @param grid_pol object of class sf, sfc, sfg or SpatialPolygons. Grid hexagones or squares.The shapefile must be in a projected coordinate system.
#' @param grid_id character. Column name of the grid ID.
#' @param grid_type character. If grid is null you can make a regular grid of "hexagonal" or "square".
#' @param cellsize numeric. Grid area (m²).
#' @param grid_boundary logical.If TRUE, the Incomplete "hexagons" or "squares" in the boundaries of the region will be discarded
#' @param clip logical. If TRUE, the new grid will be clipped to the region area. The operation time will be longer the greater the number of vertices in the polygon of the region, if it is a region with many vertices use the argument "tolerance".
#' @param tolerance numeric. If "clip" is equal to TRUE reduces the number of vertices in the region polygon.
#' @param distance list.See distancefile(). E.g.: list(type= "centroid", resistance = NULL).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 (default) that is 50 percentage of probability connection. Use in case of selecting the "PC" metric or "ProtConn".
#' @param distance_threshold numeric. Distance thresholds to establish connections.
#' @param ProtConn list. ProtConn parameters: list(attribute = "Intersected area", transboundary = NULL)
#' @param PC_IIC list. PC or IIC metric parameters: list(attribute  = NULL, attribute_weighted_area = FALSE, metric = c("IIC", "PC")).
#' @param intern logical. Show the progress of the process, default = TRUE.
#' @param parallel logical. Parallelize the function using furrr package and multiprocess plan, default = FALSE.
#' @references Matt Strimas-Mackey. http://strimas.com/spatial/hexagonal-grids/.\cr
#' Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.\url{https://doi.org/10.1016/j.ecolind.2016.12.047}.\cr
#' Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. 2006. Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology 21 (7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.
#' @examples ruta <- system.file("extdata", "WDPA_May2019_MEX-shapefile-polygons.shp", package = "Makurhini")
#' Protected_areas <- sf::read_sf(ruta)
#' ruta <- system.file("extdata",  "Region_test.shp", package = "Makurhini")
#' region <- sf::read_sf(ruta)
#' hexagons_priority <- Connect_grid(nodes = Protected_areas, region = region, grid_type = "hexagonal", cellsize = 300000000,#100km2
#'                                   distance = list(type = "centroid"), distance_thresholds = 10000, clip = TRUE,
#'                                   PC_IIC = list(metric = "PC"), intern = TRUE, parallel = TRUE)
#' hexagons_priority
#'
#' hexagons_priority_2 <- Connect_grid(nodes = Protected_areas, region = region, grid_type = "hexagonal", cellsize = 300000000,#100km2
#' distance = list(type = "centroid"), distance_thresholds = 5000, clip = TRUE, ProtConn = list(attribute ="Intersected area", transboundary = 10000),
#' intern = TRUE, parallel = TRUE)
#' hexagons_priority_2
#' @export
#' @importFrom magrittr %>%
Connect_grid <- function(nodes, region = NULL, grid_pol = NULL, grid_id = NULL,
                         grid_type = c("hexagonal", "square"),
                         cellsize = NULL, grid_boundary = FALSE,
                         clip = FALSE, tolerance = NULL,
                         distance = list(type = "centroid"),
                         probability = 0.5,
                         distance_thresholds = NULL,
                         ProtConn = list(attribute = "Intersected area",
                                         transboundary = NULL),
                         PC_IIC = list(attribute  = NULL,
                                       attribute_weighted_area = FALSE,
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
        stop("error missing cellsize(m²)")
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
    if (is.null(probability) | !is.numeric(probability)) {
      stop("error missing probability")
    }
  } else {
    if (!PC_IIC$metric %in% c("PC", "IIC")) {
      stop("Type must be either 'PC', or 'IIC'")
    }
    if (PC_IIC$metric == "PC") {
      if (is.null(probability) | !is.numeric(probability)) {
        stop("error missing probability")
      }
    }
  }

  options(warn = -1)
  ttt.2 <- getwd()
  make_grid <- function(x, type, cell_width, cell_area, clip = FALSE, tolerance = tolerance, grid_boundary = grid_boundary) {
    if (class(x)[1] == "sf") {
      x <- sf::st_zm(x)
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
    ext <- as(raster::extent(x) + cell_width, "SpatialPolygons")
    raster::projection(ext) <- raster::projection(x)

    if (type == "square") {
      g <- raster::raster(ext, resolution = cell_width)
      g <- as(g, "SpatialPolygons") %>% sf::st_as_sf() %>%
        sf::st_cast("POLYGON")
    } else if (type == "hexagonal") {
      g <- sp::spsample(ext, type = "hexagonal", cellsize = cell_width,
                    offset = c(0, 0))
      g <- sp::HexPoints2SpatialPolygons(g, dx = cell_width) %>%
        sf::st_as_sf() %>% sf::st_cast("POLYGON")
    } else {
      stop("Select a grid_type = c('square', 'hexagonal')")
    }

    x2 <- sf::st_as_sf(x) %>% sf::st_cast("POLYGON")
    if (mapview::npts(x2) > 100) {
      region_1 <- tryCatch(sf::st_simplify(x2, dTolerance = ((raster::extent(x)[4]-raster::extent(x)[3])*2)/100,
                                       preserveTopology = TRUE), error = function(err) err)
      if (inherits(region_1, "error")) {
        x2 <- sf::st_buffer(x2, dist = 0)
        x2 <- sf::st_simplify(x2, dTolerance = ((raster::extent(x)[4]-raster::extent(x)[3])*2)/100, preserveTopology = TRUE)
      } else {
        x2 <- region_1
      }
    }
    x2 <- tryCatch(x2 %>% sf::st_cast("POLYGON") %>% rmapshaper::ms_dissolve() %>%
                     sf::st_buffer(dist = 0), error = function(err) err)
    if (inherits(x2, "error")) {
      x2 <- tryCatch(x2 %>% sf::st_buffer(., dist = 0) %>%
                       rmapshaper::ms_dissolve() %>% sf::st_buffer(., dist = 0), error = function(err) err)
    }
####
    if (isTRUE(clip)) {
      x2 <- sf::st_as_sf(x) %>% sf::st_cast("POLYGON")
      if(is.null(tolerance)){
        if (mapview::npts(x2) > 100) {
          region_1 <- tryCatch(sf::st_simplify(x2, dTolerance = ((raster::extent(x)[4]-raster::extent(x)[3])*0.1)/100,
                                           preserveTopology = TRUE), error = function(err) err)
            if (inherits(region_1, "error")) {
            x2 <- sf::st_buffer(x2, dist = 0)
            x2 <- sf::st_simplify(x2, dTolerance =((raster::extent(x)[4]-raster::extent(x)[3])*0.1)/10, preserveTopology = TRUE)
          } else {
            x2 <- region_1
          }
        }
      } else {
        region_1 <- tryCatch(sf::st_simplify(x2, dTolerance = tolerance,
                                         preserveTopology = TRUE), error = function(err) err)
        if (inherits(region_1, "error")) {
          x2 <- sf::st_buffer(x2, dist = 0)
          x2 <- sf::st_simplify(x2, dTolerance = tolerance, preserveTopology = TRUE)
        } else {
          x2 <- region_1
        }
      }
      g2 <- tryCatch(sf::st_intersection(g, y = x2) %>% sf::st_cast("POLYGON"),
                     error = function(err) err)

      if (inherits(g2, "error")) {
        g2 <- tryCatch(sf::st_intersection(g, y = sf::st_buffer(x2,dist = 0)) %>%
                         sf::st_cast("POLYGON"), error = function(err) err)
      }
    } else {
      g2 <- SelectbyLoc(g, x2)

      if(isTRUE(grid_boundary)){
        if (mapview::npts(x2) > 100) {
          region_1 <- tryCatch(sf::st_simplify(x2, dTolerance = ((raster::extent(x)[4]-raster::extent(x)[3])*2)/100,
                                           preserveTopology = TRUE), error = function(err) err)
          if (inherits(region_1, "error")) {
            x2 <- sf::st_buffer(x2, dist = 0)
            x2 <- sf::st_simplify(x2, dTolerance = ((raster::extent(x)[4]-raster::extent(x)[3])*2)/100, preserveTopology = TRUE)
          } else {
            x2 <- region_1
          }
        }

        g3 <- tryCatch(sf::st_intersection(g2, y = x2) %>% sf::st_cast("POLYGON"),
                       error = function(err) err)

        if (inherits(g2, "error")) {
          g3 <- tryCatch(sf::st_intersection(g2, y = sf::st_buffer(x2,dist = 0)) %>%
                           sf::st_cast("POLYGON"), error = function(err) err)
        }

        g3$area <- as.numeric(sf::st_area(g3))
        g3 <- plot(g3[which(g3$area == cell_area), ])
      }
    }
    row.names(g2) <- as.character(1:nrow(g2))
    return(g2)
  }

  ###
  if (is.null(grid_pol) & !is.null(region)) {
    x_grid <- make_grid(x = region, type = grid_type, cell_area = cellsize, clip = clip, tolerance = tolerance, grid_boundary = grid_boundary)
  } else {
    if (class(grid_pol)[1] == "SpatialPolygonsDataFrame") {
      x_grid <- sf::st_as_sf(grid_pol) %>% sf::st_cast("POLYGON")
    } else {
      x_grid <- grid_pol
    }
  }

  ###
  x_grid$TempID <- 1:nrow(x_grid)
  x_grid <- x_grid["TempID"]
  x_grid <- base::split(x_grid, x_grid$TempID)
  if (missing(PC_IIC) & !missing(ProtConn)) {
    if (isTRUE(parallel)) {
      future::plan(strategy = future::multiprocess)
      resultado_1 <- tryCatch(furrr::future_map(x_grid, function(x) {
        test1 <- ProtConnCLa(nodes = nodes, region = x,
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
      pb <- dplyr::progress_estimated(length(x_grid), 0)
      resultado_1 <- tryCatch(purrr::map(x_grid, function(x) {
        if (isTRUE(intern)) {
          pb$tick()$print()
        }
        test1 <- ProtConnCLa(nodes = nodes, region = x,
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
      nodes1 <- sf::st_cast(nodes1, "POLYGON")
      future::plan(strategy = future::multiprocess)
      resultado_1 <- tryCatch(furrr::future_map(x_grid, function(x) {
        if (is.null(PC_IIC$attribute)) {
          nodes2 <- SelectbyLoc(nodes1, x, id = NULL,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]
          if (nrow(nodes2) > 1) {
            nodes3 <- sf::st_intersection(nodes2, x) %>%
              rmapshaper::ms_dissolve(.) %>% sf::st_buffer(., dist = 0)
            nodes2_r <- nodes2
            nodes2 <- sf::st_cast(nodes3, "POLYGON") %>%
              sf::st_buffer(., dist = 0)
            nodes2$attrib <- as.numeric(sf::st_area(nodes2)) * 1e-04

            if (nrow(nodes2) == 1) {
              nodes2 <- (nodes2$attrib * 100)/(as.numeric(sf::st_area(x)) * 1e-04)
            } else if (nrow(nodes2) == 0) {
              nodes2 = NULL
            }

          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- (nodes2$Area2 * 100)/(as.numeric(sf::st_area(x)) * 1e-04)
          } else {
            nodes2 = NULL
          }
        } else {
          nodes2 <- SelectbyLoc(nodes1, x, id = NULL,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]

          if (nrow(nodes2) > 1) {
            nodes3 <- sf::st_intersection(nodes2, x) %>%
              rmapshaper::ms_dissolve(., field = PC_IIC$attribute) %>%
              sf::st_buffer(., dist = 0)

            nodes3 <- sf::st_cast(nodes3, "POLYGON") %>%
              sf::st_buffer(., dist = 0)

            if (nrow(nodes3) > 1) {
              nodes3$IDTemp <- 1:nrow(nodes3)
              d <- distancefile(nodes3, id = "IDTemp", type = "edge")

              if (unique(d$Distance) == 0) {
                nodes2_r <- nodes2
                nodes2 <-(as.numeric(sf::st_area(nodes3)) * 1e-04)
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes2 <- nodes3[ PC_IIC$attribute][[1]] * nodes2
                }
                nodes2 <- sum((nodes2 * 100)/(as.numeric(sf::st_area(x)) * 1e-04))
              } else {
                nodes3$attrib <- nodes3[PC_IIC$attribute][[1]]
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes3$attrib <- nodes3$attrib * (as.numeric(sf::st_area(nodes3)) * 1e-04)
                }
                nodes2 <- nodes3
              }
            } else {
              nodes2_r <- nodes2
              nodes2 <-(as.numeric(sf::st_area(nodes3)) * 1e-04)
              if (isTRUE(PC_IIC$attribute_weighted_area)) {
                nodes2 <- nodes3[ PC_IIC$attribute][[1]] * nodes2
              }
              nodes2 <- sum((nodes2 * 100)/(as.numeric(sf::st_area(x)) * 1e-04))
            }
          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- nodes2$PercPi
            nodes3 <- nodes2[PC_IIC$attribute][[1]]
            if (isTRUE(PC_IIC$attribute_weighted_area)) {
              nodes3 <- nodes3 * nodes2$Area2
            }
            nodes2 <- (nodes3 * 100)/(as.numeric(sf::st_area(x)) * 1e-04)
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
            x3 <- c(NA, sum(as.numeric(sf::st_area(nodes2_r))) * 1e-04, 100) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(PC_IIC$metric, "EC", "Normalized_EC")
            x3$PArea <- if(((sum(nodes2_r$Area2 * 100))/(as.numeric(sf::st_area(x)) * 1e-04)) > 100){
              100
            } else {
              sum((nodes2_r$Area2 * 100))/(as.numeric(sf::st_area(x)) * 1e-04)
            }

            x3$EC <- if (sum(nodes2_r$Area2) > (as.numeric(sf::st_area(x)) * 1e-04)) {
              (as.numeric(sf::st_area(x)) * 1e-04)
            } else {
              sum(nodes2_r$Area2)
            }
            x3$Normalized_EC <- x3$PArea
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          }
        } else {
          test <- dConnectivity(nodes = nodes2, id = NULL,
                                attribute = "attrib", restauration = NULL,
                                distance = distance, metric = PC_IIC$metric,
                                probability = probability, distance_thresholds = distance_thresholds,
                                overall = TRUE, LA = as.numeric(sf::st_area(x)) *
                                  1e-04, write = NULL)
          test <- test[[1]][[2]]
          x2 <- as.data.frame(test[, 2])
          x3 <- as.data.frame(t(x2))
          x3 <- x3[, c(2:3)]
          names(x3) <- c("EC", PC_IIC$metric)
          x3$Normalized_EC <- (x3$EC * 100)/(as.numeric(sf::st_area(x)) *
                                               1e-04)
          x3$PArea <- (sum(nodes2$attrib) * 100)/(as.numeric(sf::st_area(x)) *
                                                    1e-04)
          x3$TempID <- paste0(x[["TempID"]])
          x3 <- x3[, c(4, 2, 1, 3, 5)]
        }
        return(x3)
      }, .progress = intern), error = function(err) err)
    } else {
      pb <- dplyr::progress_estimated(length(x_grid), 0)
      nodes1 <- nodes
      nodes1 <- sf::st_cast(nodes1, "POLYGON")
      resultado_1 <- tryCatch(purrr::map(x_grid, function(x) {
        if (isTRUE(intern)) {
          pb$tick()$print()
        }
        if (is.null(PC_IIC$attribute)) {
          nodes2 <- SelectbyLoc(nodes1, x, id = NULL,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]
          if (nrow(nodes2) > 1) {
            nodes3 <- sf::st_intersection(nodes2, x) %>%
              rmapshaper::ms_dissolve(.) %>% sf::st_buffer(., dist = 0)
            nodes2_r <- nodes2
            nodes2 <- sf::st_cast(nodes3, "POLYGON") %>%
              sf::st_buffer(., dist = 0)
            nodes2$attrib <- as.numeric(sf::st_area(nodes2)) * 1e-04

            if (nrow(nodes2) == 1) {
              nodes2 <- (nodes2$attrib * 100)/(as.numeric(sf::st_area(x)) * 1e-04)
            } else if (nrow(nodes2) == 0) {
              nodes2 = NULL
            }

          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- (nodes2$Area2 * 100)/(as.numeric(sf::st_area(x)) * 1e-04)
          } else {
            nodes2 = NULL
          }
        } else {
          nodes2 <- SelectbyLoc(nodes1, x, id = NULL,
                                selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]

          if (nrow(nodes2) > 1) {
            nodes3 <- sf::st_intersection(nodes2, x) %>%
              rmapshaper::ms_dissolve(., field = PC_IIC$attribute) %>%
              sf::st_buffer(., dist = 0)

            nodes3 <- sf::st_cast(nodes3, "POLYGON") %>%
              sf::st_buffer(., dist = 0)

            if (nrow(nodes3) > 1) {
              nodes3$IDTemp <- 1:nrow(nodes3)
              d <- distancefile(nodes3, id = "IDTemp", type = "edge")

              if (unique(d$Distance) == 0) {
                nodes2_r <- nodes2
                nodes2 <-(as.numeric(sf::st_area(nodes3)) * 1e-04)
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes2 <- nodes3[ PC_IIC$attribute][[1]] * nodes2
                }
                nodes2 <- sum((nodes2 * 100)/(as.numeric(sf::st_area(x)) * 1e-04))
              } else {
                nodes3$attrib <- nodes3[PC_IIC$attribute][[1]]
                if (isTRUE(PC_IIC$attribute_weighted_area)) {
                  nodes3$attrib <- nodes3$attrib * (as.numeric(sf::st_area(nodes3)) * 1e-04)
                }
                nodes2 <- nodes3
              }
            } else {
              nodes2_r <- nodes2
              nodes2 <-(as.numeric(sf::st_area(nodes3)) * 1e-04)
              if (isTRUE(PC_IIC$attribute_weighted_area)) {
                nodes2 <- nodes3[PC_IIC$attribute][[1]] * nodes2
              }
              nodes2 <- sum((nodes2 * 100)/(as.numeric(sf::st_area(x)) * 1e-04))
            }
          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- nodes2$PercPi
            nodes3 <- nodes2[PC_IIC$attribute][[1]]
            if (isTRUE(PC_IIC$attribute_weighted_area)) {
              nodes3 <- nodes3 * nodes2$Area2
            }
            nodes2 <- (nodes3 * 100)/(as.numeric(sf::st_area(x)) * 1e-04)
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
            x3 <- c(NA, sum(as.numeric(sf::st_area(nodes2_r))) * 1e-04, 100) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(PC_IIC$metric, "EC", "Normalized_EC")
            x3$PArea <- if(((sum(nodes2_r$Area2 * 100))/(as.numeric(sf::st_area(x)) * 1e-04)) > 100){
              100
            } else {
              sum((nodes2_r$Area2 * 100))/(as.numeric(sf::st_area(x)) * 1e-04)
            }

            x3$EC <- if (sum(nodes2_r$Area2) > (as.numeric(sf::st_area(x)) * 1e-04)) {
              (as.numeric(sf::st_area(x)) * 1e-04)
            } else {
              sum(nodes2_r$Area2)
            }
            x3$Normalized_EC <- x3$PArea
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          }
        } else {
          test <- dConnectivity(nodes = nodes2, id = NULL,
                                attribute = "attrib", restauration = NULL,
                                distance = distance, metric = PC_IIC$metric,
                                probability = probability, distance_thresholds = distance_thresholds,
                                overall = TRUE, LA = as.numeric(sf::st_area(x)) *
                                  1e-04, write = NULL)
          test <- test[[1]][[2]]
          x2 <- as.data.frame(test[, 2])
          x3 <- as.data.frame(t(x2))
          x3 <- x3[, c(2:3)]
          names(x3) <- c("EC", PC_IIC$metric)
          x3$Normalized_EC <- (x3$EC * 100)/(as.numeric(sf::st_area(x)) *
                                               1e-04)
          x3$PArea <- (sum(nodes2$attrib) * 100)/(as.numeric(sf::st_area(x)) *
                                                    1e-04)
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
    x_grid2 <- base::merge(x_grid2, result_2, by = "TempID")
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
