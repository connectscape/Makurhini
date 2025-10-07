#' Estimate the ECA, dA and dECA in a regular grid
#'
#' Use the function to compute the ECA, dECA, dA indexes in a regular grid.
#' @param nodes \code{list}. List of objects containing nodes (e.g., habitat patches or fragments) of each time to analyze information. Nodes are spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.
#' @param nodes_names \code{character}. (\emph{optional, default =} \code{NULL}). Name of each time or scenario used in the nodes parameter
#' @param attribute \code{character}. If \code{NULL} the area of the nodes will be used as the node attribute. The unit of area can be selected using the \code{area_unit} parameter.
#'  Specify \bold{the name of the column} containing the attribute for the nodes. The column name must be present in each element of the node list.\cr
#' @param area_unit \code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.
#' @param region object of class \code{sf}, \code{SpatialPolygonsDataFrame}. Polygon delimiting the region or study area. It must be
#'  in a projected coordinate system.
#' @param grid \code{list} or object of class \code{sf}, \code{SpatialPolygonsDataFrame}.
#' Use this parameter to generate a grid indicating its characteristics in a \code{list} (see \link[Makurhini]{get_grid}) or enter the name of an sf class \code{sf} or \code{SpatialPolygonsDataFrame} with the grid whose coordinate system must be the same as that of the \code{nodes}.
#'  Example for generating 100 km\out{<sup>2</sup>} hexagons:\cr
#' \code{list(hexagonal = TRUE, cellsize = unit_convert(100, "km2", "m2"), grid_boundary = FALSE, clip = FALSE, tolerance = NULL)}.
#' @param distance  A \code{list} to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. It must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. To see more arguments see the \link[Makurhini]{distancefile} function.
#'  You can place a \code{list} with resistances where there must be one resistance for each time/scenario of patches in the \code{nodes} parameter, for example,
#'  if nodes has a list with two time patches then you can use two resistances one for time 1 and one for time 2:\cr
#'  \code{distance(type = "least-cost", resistance = list(resistanceT1, resistanceT2))}.
#' @param metric A \code{character} indicating the connectivity metric to use: \code{"PC"} (the default and recommended) to calculate the probability of connectivity index, and \code{"IIC"} to calculate the binary integral index of connectivity.
#' @param probability A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the \code{"PC"} metric. If \code{probability = NULL}, then a probability of 0.5 will be used.
#' @param distance_threshold A \code{numeric} indicating the dispersal distance (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range. Can be the same length as the \code{distance_thresholds} parameter.
#' @param threshold \code{numeric}. Pairs of nodes with a distance value greater than this threshold will be discarded in the analysis which can speed up processing.
#' @param parallel \code{numeric}. Specify the number of cores to use for parallel processing, \code{default = NULL}. Parallelize the function using furrr package.
#' @param intern \code{logical}. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @return List with sf objects corresponding to the hexagons and each transition between scenarios or node times, each of the following fields:\cr\cr
#' -  Time: name of the time periods, name of the model or scenario (are taken from the name of the elements of the list of nodes or the plot argument)\cr
#' -  ECA.i: Equivalent Connected Area or Equivalent Connectivity for time i (first time/scenery in the comparison)\cr
#' -  ECA.j: Equivalent Connected Area or Equivalent Connectivity for time j (second time/scenery in the comparison)\cr
#' -  ECA.i: Normalized_ECA (% of LA) or relative connectivity for time i\cr
#' -  ECA.j: Normalized_ECA (% of LA) or relative connectivity for time j\cr
#' -  dA: delta Area between times (percentage)\cr
#' -  dECA: delta ECA between times (percentage)\cr
#' -  rECA: relativized ECA (dECA/dA). According to Liang et al. (2021) "an rECA value greater than 1 indicates that habitat changes result in a
#' disproportionately large change in habitat connectivity, while a value lower than 1 indicates connectivity
#' changes due to random habitat changes (Saura et al. 2011; Dilts et al. 2016)".\cr
#' -  dA/dECA comparisons: comparisons between dA and dECA\cr
#' -  Type of change: Type of change using the dECAfun() and the difference between dA and dECA.\cr
#' @references \url{www.conefor.org}\cr\cr
#' - Saura, S., Estreguil, C., Mouton, C., & Rodríguez-Freire, M. (2011). Network analysis to assess landscape connectivity trends: Application to European forests (1990-2000). Ecological Indicators, 11(2), 407–416.
#' https://doi.org/10.1016/j.ecolind.2010.06.011 \cr
#'  Herrera, L. P., Sabatino, M. C., Jaimes, F. R., & Saura, S. (2017). Landscape connectivity and the role of small habitat patches as stepping stones: an assessment of the grassland biome in South America. Biodiversity and Conservation, 26(14), 3465–3479.
#' https://doi.org/10.1007/s10531-017-1416-7\cr
#' - Liang, J., Ding, Z., Jiang, Z., Yang, X., Xiao, R., Singh, P. B., ... & Hu, H. (2021). Climate change, habitat connectivity, and conservation gaps: a case study of four ungulate species endemic to the Tibetan Plateau. Landscape Ecology, 36(4), 1071-1087.\cr
#' - Dilts TE, Weisberg PJ, Leitner P, Matocq MD, Inman RD, Nussear KE, Esque TC (2016) Multi-scale connectivity and graph theory highlight critical areas for conservation under climate change. Ecol Appl 26:1223–1237
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(sf)
#'
#' # Four times (T1.2, T2.3, T3.4)
#' data("list_forest_patches", package = "Makurhini")
#' data("study_area", package = "Makurhini")
#' class(list_forest_patches)
#'
#' hexagons_dECA <- MK_dECA_grid(nodes = list_forest_patches,
#'                               nodes_names = c("T1", "T2", "T3", "T4"),
#'                               region = study_area,
#'                               area_unit = "ha",
#'                               metric = "IIC",
#'                               grid = list(hexagonal = TRUE,
#'                                                 cellsize = unit_convert(100, "km2", "m2")),
#'                               distance_threshold = 3000,
#'                               probability = 0.5,
#'                               distance = list(type = "centroid"),
#'                               parallel = 4,
#'                               intern = TRUE)
#' names(hexagons_dECA)#List of lenght 3, where each element is a transition.
#' plot(hexagons_dECA$result_T1.T2 ["dECA"], breaks = "quantile")
#' plot(hexagons_dECA$result_T3.T4["Type.Change"], key.pos = 1)
#'}
#' @return The function returns a list comprising elements corresponding to discrete periods. Thus, if a list of nodes contains three scenarios or times, the function returns a list with two elements. The first element corresponds to the transition between scenarios 1 and 2, and it will include the dECA value for that period. The second element of the list corresponds to the transition between scenarios 2 and 3, and it will include the dECA value for that period.
#' @importFrom sf st_as_sf st_area st_intersection st_is_empty
#' @importFrom future plan multicore multisession availableCores
#' @importFrom furrr future_map_dfr
#' @importFrom purrr map_dfr map_chr
#' @importFrom rmapshaper ms_explode
#' @export

MK_dECA_grid <- function(nodes,
                         nodes_names = NULL,
                         attribute = NULL,
                         region = NULL,
                         grid = list(hexagonal = TRUE,
                                     cellsize = NULL, grid_boundary = FALSE,
                                     clip = FALSE, tolerance = NULL),
                         area_unit = "m2",
                         distance = list(type = "centroid"),
                         metric = "IIC",
                         distance_threshold = NULL,
                         threshold,
                         probability = NULL,
                         parallel = NULL,
                         intern = TRUE){
  options(warn = -1)

  if(isFALSE(parallel)){
    parallel <- NULL
  }

  if(isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
  }

  if(!is.null(parallel)){
    if(!is.numeric(parallel)){
      stop("if you use parallel argument then you need a numeric value")
    }
  }

  message("Step 1. Reviewing parameters")
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (!metric %in% c("IIC", "PC")) {
    stop("Type must be either 'IIC', or 'PC'")
  }

  if (isTRUE(unique(metric == c("IIC", "PC")))) {
    metric = "IIC"
  }

  if (metric == "PC") {
    if (!is.null(probability) & !is.numeric(probability)) {
      stop("error missing probability")
    }
  }

  if (is.null(distance_threshold)) {
    stop("error missing numeric distance threshold(s)")
  }

  if(class(grid)[1] == "list"){
    message("Step 2. Grid processing")
    base_grid <- get_grid(region = region,
                          hexagonal = grid$hexagonal,
                          cellsize = grid$cellsize,
                          grid_boundary = grid$grid_boundary,
                          clip = grid$clip,
                          tolerance = grid$tolerance)
  } else {
    if(any(class(grid)[1] == "sf" | class(grid)[1] == "SpatialPolygonsDataFrame")){
      base_grid <- get_grid(grid_pol = grid)
    } else {
      stop("You need the grid parameter")
    }
  }

  if(is.null(nodes_names)){
    nodes_names <- paste0("T", 1:length(nodes))
  }

  if(class(base_grid)[1] != "grid"){
    stop("error making the grid")
  }

  if (isTRUE(intern)) {
    message("Step 3. Processing metric on the grid. Progress estimated:")
  } else {
    message("Step 3. Processing metric on the grid")
  }

  if (!is.null(parallel)){
    works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}
    if(.Platform$OS.type == "unix") {
      strat <- future::multicore
    } else {
      strat <- future::multisession
    }
    plan(strategy = strat, gc = TRUE, workers = works)
    nodes <- future_map(nodes, function(x){
      if(class(x)[[1]] == "SpatialPolygonsDataFrame"){
        x <- st_as_sf(x)
      }
      return(x)
    })
    result_1 <- future_map_dfr(1:nrow(base_grid@grid), function(x){
      x.1 <- base_grid@grid[x,]
      dECA.1 <- data.frame("Time" = paste0(nodes_names[1:(length(nodes)-1)], ".", nodes_names[2:length(nodes)]),
                           "ECA.i" = 0,
                           "ECA.j" = 0,
                           "ECAnorm.i" = 0,
                           "ECAnorm.j" = 0,
                           "dA" = 0,
                           "dECA" = 0,
                           "rECA" = 0,
                           "comparisons" = "NA",
                           "Type.Change" = "NA")

      LA <- unit_convert(st_area(x.1), "m2", area_unit)
      nodes.1 <- lapply(nodes, function(i){
        i.1 <- suppressWarnings(st_intersection(i, x.1)); i.1 <- i.1[which(!st_is_empty(i.1)),]

        if(nrow(i.1) > 0){
          i.1 <- ms_explode(i.1); i.1$IdTemp <- 1:nrow(i.1); i.1 <- i.1["IdTemp"]
          i.1 <- i.1[which(!st_is_empty(i.1)),]
        } else {
          i.1 <- "NA"
        }
        return(i.1)
      })

      if(!is.null(distance$resistance)){
        distance2 <- distance
        if(class(distance2$resistance)[1] == "list"){
          x.0 <- st_buffer(x.1, res(distance2$resistance[[1]])[1]*10)
          distance2$resistance <- raster::crop(distance2$resistance, x.0)
        } else {
          distance2$resistance <- tryCatch(lapply(distance2$resistance, function(x){
            x.0 <- st_buffer(x.1, res(distance2$resistance[[1]])[1]*10)
            x.1 <- raster::crop(distance2$resistance, x.0)
            return(x.1)
          }), error = function(err)err)
        }
      } else {
        distance2 <- distance
      }
      dECA.2 <- map_dfr(1:(length(nodes.1)-1), function(i){
        i.1 <- dECA.1[i,]; i.1$t <- i.1$Time; i.1 <- i.1[,2:ncol(i.1)]

        if(!is.character(nodes.1[[i]])){
          if(is.character(nodes.1[[i+1]])){ #shp-NA
            if(nrow(nodes.1[[i]]) > 2){
              if(!is.null(distance2$resistance)){
                if(class(distance2$resistance)[1] == "list"){
                  distance2$resistance <- distance2$resistance[[i]]
                }
              }

              ECAi <- MK_dPCIIC(nodes = nodes.1[[i]],
                                attribute = attribute,
                                area_unit = area_unit,
                                distance = distance2,
                                metric = metric,
                                probability = probability,
                                distance_thresholds = distance_threshold,
                                onlyoverall = TRUE,
                                LA = LA, intern = FALSE)
              ECAi <- ECAi[2,2]
            } else {
              ECAi <- sum(unit_convert(st_area(nodes.1[[i]]), "m2", area_unit))
            }
            i.1[,1] <- ECAi; i.1[,5:6] <- -100; i.1[,3] <- (i.1[,1]*100)/LA
            i.1[,7] <- 1; i.1[,8] <- dECAfun(-100, -100); i.1[,9] <- dECAfun2(-100, -100)
          } else {  #shp-shp
            dECA.2 <- tryCatch(MK_dECA(nodes = nodes.1[i:(i+1)],
                                       attribute = attribute,
                                       area_unit = area_unit,
                                       distance = distance2,
                                       metric = metric,
                                       probability = probability,
                                       distance_thresholds = distance_threshold,
                                       LA = LA,
                                       plot = FALSE, parallel = NULL,
                                       write = NULL, intern = FALSE), error = function(err)err)
            i.1[,5:9] <- dECA.2[2,8:12]; i.1[1:2] <- dECA.2[[5]]
            i.1[1:2] <- dECA.2[[5]]; i.1[3:4] <- dECA.2[[6]]
          }

        } else {
          if(is.character(nodes.1[[i+1]])){
            i.1[,1:7] <- NA
          } else {
            i.1[,2] <- LA; i.1[,4] <- 100
            i.1[,5:6] <- 100; i.1[,7] <- 1; i.1[,8] <- dECAfun(100, 100); i.1[,9] <- dECAfun2(100, 100)
          }
        }
        names(i.1)[ncol(i.1)] <- "Time"; i.1 <- i.1[,c(ncol(i.1), 1:(ncol(i.1)-1))]
        return(i.1)
      })
      return(dECA.2)
    }, .progress = intern)
    close_multiprocess(works)
  } else {
    nodes <- lapply(nodes, function(x){
      if(class(x)[[1]] == "SpatialPolygonsDataFrame"){
        x <- st_as_sf(x)
      }
      return(x)
    })
    result_1 <- map_dfr(1:nrow(base_grid@grid), function(x){
      x.1 <- base_grid@grid[x,]
      dECA.1 <- data.frame("Time" = paste0(nodes_names[1:(length(nodes)-1)], ".", nodes_names[2:length(nodes)]),
                           "ECA.i" = 0,
                           "ECA.j" = 0,
                           "ECAnorm.i" = 0,
                           "ECAnorm.j" = 0,
                           "dA" = 0,
                           "dECA" = 0,
                           "rECA" = 0,
                           "comparisons" = "NA",
                           "Type.Change" = "NA")

      LA <- unit_convert(st_area(x.1), "m2", area_unit)
      nodes.1 <- lapply(nodes, function(i){
        i.1 <- suppressWarnings(st_intersection(i, x.1)); i.1 <- i.1[which(!st_is_empty(i.1)),]

        if(nrow(i.1) > 0){
          i.1 <- ms_explode(i.1); i.1$IdTemp <- 1:nrow(i.1); i.1 <- i.1["IdTemp"]
          i.1 <- i.1[which(!st_is_empty(i.1)),]
        } else {
          i.1 <- "NA"
        }
        return(i.1)
      })

      if(!is.null(distance$resistance)){
        distance2 <- distance
        if(class(distance2$resistance)[1] == "list"){
          x.0 <- st_buffer(x.1, res(distance2$resistance[[1]])[1]*10)
          distance2$resistance <- raster::crop(distance2$resistance, x.0)
        } else {
          distance2$resistance <- tryCatch(lapply(distance2$resistance, function(x){
            x.0 <- st_buffer(x.1, res(distance2$resistance[[1]])[1]*10)
            x.1 <- raster::crop(distance2$resistance, x.0)
            return(x.1)
          }), error = function(err)err)
        }
      } else {
        distance2 <- distance
      }
      dECA.2 <- map_dfr(1:(length(nodes.1)-1), function(i){
        i.1 <- dECA.1[i,]; i.1$t <- i.1$Time; i.1 <- i.1[,2:ncol(i.1)]

        if(!is.character(nodes.1[[i]])){
          if(is.character(nodes.1[[i+1]])){ #shp-NA
            if(nrow(nodes.1[[i]]) > 2){
              if(!is.null(distance2$resistance)){
                if(class(distance2$resistance)[1] == "list"){
                  distance2$resistance <- distance2$resistance[[i]]
                }
              }

              ECAi <- MK_dPCIIC(nodes = nodes.1[[i]],
                                attribute = attribute,
                                area_unit = area_unit,
                                distance = distance2,
                                metric = metric,
                                probability = probability,
                                distance_thresholds = distance_threshold,
                                onlyoverall = TRUE,
                                LA = LA, intern = FALSE)
              ECAi <- ECAi[2,2]
            } else {
              ECAi <- sum(unit_convert(st_area(nodes.1[[i]]), "m2", area_unit))
            }
            i.1[,1] <- ECAi; i.1[,5:6] <- -100; i.1[,3] <- (i.1[,1]*100)/LA
            i.1[,7] <- 1; i.1[,8] <- dECAfun(-100, -100); i.1[,9] <- dECAfun2(-100, -100)
          } else {  #shp-shp
            dECA.2 <- tryCatch(MK_dECA(nodes = nodes.1[i:(i+1)],
                                       attribute = attribute,
                                       area_unit = area_unit,
                                       distance = distance2,
                                       metric = metric,
                                       probability = probability,
                                       distance_thresholds = distance_threshold,
                                       LA = LA,
                                       plot = FALSE, parallel = NULL,
                                       write = NULL, intern = FALSE), error = function(err)err)
            i.1[,5:9] <- dECA.2[2,8:12]; i.1[1:2] <- dECA.2[[5]]
            i.1[1:2] <- dECA.2[[5]]; i.1[3:4] <- dECA.2[[6]]
          }

        } else {
          if(is.character(nodes.1[[i+1]])){
            i.1[,1:7] <- NA
          } else {
            i.1[,2] <- LA; i.1[,4] <- 100
            i.1[,5:6] <- 100; i.1[,7] <- 1; i.1[,8] <- dECAfun(100, 100); i.1[,9] <- dECAfun2(100, 100)
          }
        }
        names(i.1)[ncol(i.1)] <- "Time"; i.1 <- i.1[,c(ncol(i.1), 1:(ncol(i.1)-1))]
        return(i.1)
      })
      return(dECA.2)
    }, .progress = intern)
  }
  result_1 <- lapply(unique(result_1$Time), function(x){
    x.1 <- result_1[result_1$Time == x,]; x.1 <- cbind(base_grid@grid, x.1)
    x.1$IdTemp  <- NULL
    return(x.1)
  })
  names(result_1) <- map_chr(result_1, function(x){paste0("result_", x[["Time"]][1])})
  if(length(result_1) == 1){
    result_1 <- result_1[[1]]
  }
  return(result_1)
}

