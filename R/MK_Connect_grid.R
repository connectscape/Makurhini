#' Connectivity indexes in a regular grid
#'
#' Use the function to compute the Protected Connected (ProtConn)or EC, PC, IIC indexes on a grid.
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
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. Nodes shapefile, the shapefile must be
#'  in a projected coordinate system.
#' @param region object of class sf, sfc, sfg or SpatialPolygons. Region shapefile, the shapefile must be
#'  in a projected coordinate system.
#' @param attribute character. Select the nodes attribute: "Intersected area" = Available if the metric argument
#' is equal to ProtConn, it corresponds to intersected Protected areas (default); or
#' another specific column name with the nodes attribute, ideally this attribute mus be an area-weighted index,
#'  otherwise the interpretation of the ProtConn or PC metric may change.
#' @param thintersect numeric. Only available if you selected ProtConn as metric. Threshold of intersection in percentage allowed to select or not a target geometry.
#'  Example, if thintersect is equal to 90 then a node will be selected only if the intersection between the node and
#'   the region is >= 90 percentage. If NULL, thintersect will be 0 (default). (see, "Makurhini::MK_ProtConn")
#' @param area_unit character. Attribute area units. You can set an area unit, "Makurhini::unit_covert()" compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to hectares "ha".
#' @param metric character. Choose a connectivity metric: "ProtConn" Protected Connected Land or "PC" Probability of conectivity considering maximum product probabilities.
#' @param distance list. See distancefile(). Example, list(type= "centroid", resistance = NULL).
#' @param distance_threshold numeric. Distance threshold to establish connections (meters).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5
#' (default) that is 50 percentage of probability connection. Use in case of selecting the "PC"
#' metric or "ProtConn".
#' @param transboundary numeric. Buffer to select polygones in a second round, their attribute value = 0,
#' see  "Makurhini::MK_ProtConn".
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
#' data("regions", package = "Makurhini")
#' ecoregion <- regions[2,]
#' plot(ecoregion, col="blue")
#'
#'hexagons_priority <- MK_Connect_grid(grid_type = "hexagonal",
#'                                     cellsize = 10000, grid_boundary = FALSE,
#'                                     clip = FALSE, nodes = Protected_areas, region = ecoregion,
#'                                     attribute = "Intersected area", thintersect = NULL,
#'                                     area_unit = "ha", metric = "ProtConn",
#'                                     distance = list(type = "centroid"),
#'                                     distance_threshold = 3000,
#'                                     probability = 0.5, transboundary = 6000,
#'                                     intern = TRUE, parallel = FALSE)
#' hexagons_priority
#' plot(hexagons_priority["ProtConn"])
#' }
#' @importFrom magrittr %>%
#' @importFrom raster extent projection projection<- raster
#' @importFrom sf st_as_sf st_zm st_cast st_buffer st_area st_intersection
#' @importFrom rmapshaper ms_dissolve
#' @importFrom future plan multiprocess availableCores
#' @importFrom furrr future_map
#' @importFrom dplyr progress_estimated
#' @importFrom purrr map

MK_Connect_grid <- function(grid_pol = NULL, grid_id = NULL,
                            grid_type = c("hexagonal", "square"),
                            cellsize = NULL, grid_boundary = FALSE,
                            clip = FALSE, tolerance = NULL,
                            nodes, region = NULL,
                            attribute = NULL,
                            thintersect = NULL,
                            area_unit = "ha",
                            metric = c("ProtConn", "PC"),
                            distance = list(type = "centroid"),
                            distance_threshold = NULL,
                            probability = NULL,
                            transboundary = NULL,
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

  if (is.null(distance_threshold)) {
    stop("error missing numeric distance threshold(s)")
  }

  if (metric == "ProtConn") {
    if (!is.null(probability) & !is.numeric(probability)) {
      stop("error missing probability")
    }
  } else {
    if (metric !=  "PC") {
      stop("Type must be'PC'")
    } else {
      if (!is.null(probability) & !is.numeric(probability)) {
        stop("error missing probability")
      }
    }
  }

  if(is.null(area_unit)){
    area_unit = "ha"
  }

  if(metric == "ProtConn" & is.null(attribute)){
    attribute = "Intersected area"
  }

  if(metric == "ProtConn" & is.null(thintersect)){
    thintersect = 0
  }

  if(metric == "ProtConn" & is.character(thintersect)){
    stop("thintersect must be NULL or numeric")
  }

  options(warn = -1)

  if (class(nodes)[1] == "SpatialPolygonsDataFrame") {
    nodes <- st_as_sf(nodes) %>% st_zm()
  } else {
    nodes <- st_zm(nodes)
  }

  nodes$IdTemp <- 1:nrow(nodes)

  if (is.null(grid_pol) & !is.null(region)) {
    x_grid <- make_grid(x = region, type = grid_type, cell_area = unit_convert(cellsize, "km2", "m2"), clip = clip, tolerance = tolerance, grid_boundary = grid_boundary)
  } else {
    if(class(grid_pol)[1] == "SpatialPolygonsDataFrame"){
      x_grid <- st_as_sf(grid_pol) %>% st_zm() %>% st_cast("POLYGON")
    } else {
      x_grid <- grid_pol %>% st_zm() %>% st_cast("POLYGON")
    }
  }

  select_distance <- max(c(transboundary, distance_threshold)) * 2
  mask.1 <- ms_dissolve(x_grid) %>% st_buffer(., dist = select_distance)

  nodes.2 <- tryCatch(over_poly(x = nodes, y = mask.1), error = function(err)err)

  if (inherits(nodes.2, "error")){
    nodes <- st_buffer(x = nodes, dist = 0)
    nodes.2 <- over_poly(x = nodes, y = mask.1)
    nodes.2 <- nodes[which(!is.na(nodes.2)),]
  } else {
    nodes.2 <- nodes[which(!is.na(nodes.2)),]
  }

  ###grid to list
  x_grid$TempID <- 1:nrow(x_grid)
  x_grid <- x_grid["TempID"]
  x_grid <- split(x_grid, x_grid$TempID)
  if (metric == "ProtConn") {
    attribute <- if(is.null(attribute)){"Intersected area"} else{attribute}
    if (isTRUE(parallel)) {
      works <- as.numeric(availableCores())-1
      plan(strategy = multiprocess, gc = TRUE, workers = works)
      resultado_1 <- tryCatch(future_map(x_grid, function(x) {
        ProtConn.1 <- MK_ProtConn(nodes = nodes.2, region = x,
                                  area_unit = area_unit,
                                  thintersect = thintersect, attribute = attribute,
                                  distance = distance, distance_thresholds = distance_threshold,
                                  probability = probability, transboundary = transboundary,
                                  LA = NULL, plot = FALSE, dPC = FALSE, write = NULL,
                                  intern = FALSE)

        n <- as.vector(ProtConn.1[[3]])
        x.2 <- as.data.frame(t(ProtConn.1[[4]]))
        names(x.2) <- n
        x.2$TempID <- paste0(x[["TempID"]])
        return(x.2)
      }, .progress = intern), error = function(err) err)
      close_multiprocess(works)
    } else {
      pb <- progress_estimated(length(x_grid), 0)
      resultado_1 <- tryCatch(map(x_grid, function(x) {
        if (isTRUE(intern)) {
          pb$tick()$print()
        }

        ProtConn.1 <- MK_ProtConn(nodes = nodes.2, region = x,
                                  area_unit = area_unit,
                                  thintersect = thintersect, attribute = attribute,
                                  distance = distance, distance_thresholds = distance_threshold,
                                  probability = probability, transboundary = transboundary,
                                  LA = NULL, plot = FALSE, dPC = FALSE, write = NULL,
                                  intern = FALSE)

        n <- as.vector(ProtConn.1[[3]])
        x.2 <- as.data.frame(t(ProtConn.1[[4]]))
        names(x.2) <- n
        x.2$TempID <- paste0(x[["TempID"]])
        return(x.2)
      }), error = function(err) err)
    }
  } else if (metric == "PC") {
    if (isTRUE(parallel)) {
      nodes.3 <- nodes.2
      nodes.3 <- st_cast(nodes.3, "POLYGON") %>% st_zm()
      nodes.3$IdTemp <- 1:nrow(nodes.3)
      works <- as.numeric(availableCores())-1
      plan(strategy = multiprocess, gc = TRUE, workers = works)
      resultado_1 <- tryCatch(future_map(x_grid, function(x) {
        area_x <- unit_convert(as.numeric(st_area(x)), "m2", area_unit)
        if (is.null(attribute)) {
          nodes2 <- MK_selectbyloc(nodes.3, x, id = NULL, area_unit = area_unit,
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
          nodes2 <- MK_selectbyloc(nodes.3, x, id = NULL, area_unit = area_unit,
                                   selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]

          if (nrow(nodes2) > 1) {
            nodes3 <- st_intersection(nodes2, x) %>%
              ms_dissolve(., field = attribute) %>%
              st_buffer(., dist = 0)

            nodes3 <- st_cast(nodes3, "POLYGON") %>%
              st_buffer(., dist = 0)

            if (nrow(nodes3) > 1) {
              nodes3$IDTemp <- 1:nrow(nodes3)
              d <- distancefile(nodes3, id = "IDTemp", type = "edge")

              if (unique(d$Distance) == 0) {
                nodes2_r <- nodes2
                nodes2 <- unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
                nodes2 <- sum((nodes2 * 100)/area_x)
              } else {
                nodes3$attrib <- nodes3[attribute][[1]]
                nodes2 <- nodes3
              }
            } else {
              nodes2_r <- nodes2
              nodes2 <-unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
              nodes2 <- sum((nodes2 * 100)/area_x)
            }
          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- nodes2$PercPi
            nodes3 <- nodes2[attribute][[1]]
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
            names(x3) <- c(metric, "EC", "Normalized_EC")
            x3$PArea <- 0
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          } else {
            x3 <- c(NA, sum(unit_convert(as.numeric(st_area(nodes2_r)), "m2", area_unit)), 100) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(metric, "EC", "Normalized_EC")

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
          test <- MK_dPCIIC(nodes = nodes2,
                            attribute = "attrib", restauration = NULL,
                            distance = distance, metric = "PC",
                            probability = probability,
                            distance_thresholds = distance_threshold,
                            overall = TRUE, LA = area_x, onlyoverall = TRUE,
                            write = NULL)
          x2 <- as.data.frame(test[,2])
          x3 <- as.data.frame(t(x2))
          x3 <- x3[, c(2:3)]
          names(x3) <- c("EC", metric)
          x3$Normalized_EC <- (x3$EC * 100)/area_x
          x3$PArea <- (sum(nodes2$attrib) * 100)/area_x
          x3$TempID <- paste0(x[["TempID"]])
          x3 <- x3[, c(4, 2, 1, 3, 5)]
        }
        return(x3)
      }, .progress = intern), error = function(err) err)
      close_multiprocess(works)

    } else {
      pb <- progress_estimated(length(x_grid), 0)
      nodes.3 <- nodes.2
      nodes.3 <- st_cast(nodes.3, "POLYGON") %>% st_zm()
      nodes.3$IdTemp <- 1:nrow(nodes.3)
      resultado_1 <- tryCatch(map(x_grid, function(x) {
        if (isTRUE(intern)) {
          pb$tick()$print()
        }
        area_x <- unit_convert(as.numeric(st_area(x)), "m2", area_unit)
        if (is.null(attribute)) {
          nodes2 <- MK_selectbyloc(nodes.3, x, id = NULL, area_unit = area_unit,
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
          nodes2 <- MK_selectbyloc(nodes.3, x, id = NULL, area_unit = area_unit,
                                   selreg = "M2", transboundary = 10)
          nodes2 <- nodes2[which(nodes2$transboundary == 1), ]

          if (nrow(nodes2) > 1) {
            nodes3 <- st_intersection(nodes2, x) %>%
              ms_dissolve(., field = attribute) %>%
              st_buffer(., dist = 0)

            nodes3 <- st_cast(nodes3, "POLYGON") %>%
              st_buffer(., dist = 0)

            if (nrow(nodes3) > 1) {
              nodes3$IDTemp <- 1:nrow(nodes3)
              d <- distancefile(nodes3, id = "IDTemp", type = "edge")

              if (unique(d$Distance) == 0) {
                nodes2_r <- nodes2
                nodes2 <- unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
                nodes2 <- sum((nodes2 * 100)/area_x)
              } else {
                nodes3$attrib <- nodes3[attribute][[1]]
                nodes2 <- nodes3
              }
            } else {
              nodes2_r <- nodes2
              nodes2 <-unit_convert(as.numeric(st_area(nodes3)), "m2", area_unit)
              nodes2 <- sum((nodes2 * 100)/area_x)
            }
          } else if (nrow(nodes2) == 1) {
            nodes2_r <- nodes2
            nodes2 <- nodes2$PercPi
            nodes3 <- nodes2[attribute][[1]]
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
            names(x3) <- c(metric, "EC", "Normalized_EC")
            x3$PArea <- 0
            x3$TempID <- paste0(x[["TempID"]])
            x3 <- x3[, c(4, 1:3, 5)]
          } else {
            x3 <- c(NA, sum(unit_convert(as.numeric(st_area(nodes2_r)), "m2", area_unit)), 100) %>% as.data.frame()
            x3 <- as.data.frame(t(x3))
            names(x3) <- c(metric, "EC", "Normalized_EC")

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
          test <- MK_dPCIIC(nodes = nodes2,
                            attribute = "attrib", restauration = NULL,
                            distance = distance, metric = "PC",
                            probability = probability,
                            distance_thresholds = distance_threshold,
                            overall = TRUE, LA = area_x, onlyoverall = TRUE,
                            write = NULL)
          x2 <- as.data.frame(test[,2])
          x3 <- as.data.frame(t(x2))
          x3 <- x3[, c(2:3)]
          names(x3) <- c("EC", metric)
          x3$Normalized_EC <- (x3$EC * 100)/area_x
          x3$PArea <- (sum(nodes2$attrib) * 100)/area_x
          x3$TempID <- paste0(x[["TempID"]])
          x3 <- x3[, c(4, 2, 1, 3, 5)]
        }
        return(x3)}), error = function(err) err)
    }
  } else {
    stop("check arguments or select one metric ProtConn or PC")
  }
  ####
  if (inherits(resultado_1, "error")) {
    stop(resultado_1)
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
  }
  return(resultado_1)
}
