#' Protected Connected (ProtConn)
#'
#' Estimate Protected Connected (ProtConn) indicator and fractions for one region.
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param region object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param thintersect numeric. Threshold of intersection in percentage allowed to select or not a target geometry.
#'  Example, if thintersect is equal to 90 then a node will be selected only if the intersection between the node and
#'   the region is >= 90 percentage. If NULL, thintersect will be 0 (default)
#' @param area_unit character. Attribute area units. You can set an area unit, "Makurhini::unit_covert()" compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to hectares "m2".
#' @param distance list. See \link[Makurhini]{distancefile}. Example, list(type= "centroid", resistance = NULL).
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param probability numeric. Probability of direct dispersal between nodes, Default, 0.5,
#'  that is 50 percentage of probability connection. If probability = NULL, then it will be the inverse of the mean dispersal distance
#' for the species (1/α; Hanski and Ovaskainen 2000).
#' @param transboundary numeric. Buffer to select polygons in a second round, their attribute value = 0, see Saura et al. 2017. You can set one transboundary value or one per each threshold distance.
#' @param transboundary_type character. Two options: "nodes" or "region". If it is "nodes" the transboundary is built from the limits of the nodes present in the region (default),
#' if "region" is selected the transboundary is built from the limits of the region.
#' @param protconn_bound logical. If TRUE then the fractions ProtUnConn[design] and ProtConn[bound] will be estimated.
#' @param LA numeric. Maximum Landscape Attribute.
#' @param geom_simplify logical. Slightly simplify the region and nodes geometries.
#' @param plot logical. Plot the main ProtConn indicators and fractions, default = FALSE.
#' @param write character. Output folder including the output file name without extension, e.g., "C:/ProtConn/Protfiles".
#' @param parallel logical. Specify the number of cores to use for parallel processing, default = NULL. Parallelize the function using furrr package and multiprocess
#' plan when there are more than ONE transboundary.
#' @param intern logical. Show the progress of the process, default = TRUE.
#' @return
#' Table with the following ProtConn values: ECA, Prot, ProtConn, ProtUnconn, RelConn, ProtUnConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within],
#'  ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot], ProtConn[Within][land], ProtConn[Contig][land],
#'  ProtConn[Unprot][land], ProtConn[Trans][land] \cr
#' \cr
#' *If plot is not NULL a list is returned with the ProtConn table and a plots.
#' @references
#' Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S., Bertzky, B., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2018). Protected area connectivity: Shortfalls in global targets and country-level priorities. Biological Conservation, 219(October 2017), 53–67.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(raster)
#' data("Protected_areas", package = "Makurhini")
#' #plot(Protected_areas, col="green")
#'
#' data("regions", package = "Makurhini")
#' region <- regions[2,]
#' plot(region, col="blue")
#'
#' test <- MK_ProtConn(nodes = Protected_areas, region = region,
#'                     area_unit = "ha",
#'                     distance = list(type= "centroid"),
#'                     distance_thresholds = c(50000, 10000),
#'                     probability = 0.5, transboundary = c(50000, 10000),
#'                     LA = NULL, plot = TRUE, parallel = NULL,
#'                     write = NULL, intern = FALSE)
#' test
#'
#' #Least-cost distances
#' HFP_Mexico <- raster(system.file("extdata", "HFP_Mexico.tif", package = "Makurhini", mustWork = TRUE))
#' mask_1 <- as(extent(Protected_areas), 'SpatialPolygons')
#' crs(mask_1) <- crs(Protected_areas)
#' mask_1 <- buffer(mask_1, 20000)
#' HFP_Mexico <- crop(HFP_Mexico, mask_1)
#'
#' #Penalizing 100 times the distance of a pixel (1 km) is too much
#' #so we will rescale the values from 1 to 10
#' HFP_Mexico <- HFP_Mexico/10
#'
#' # we have values less than 1 and that could underestimate the least cost distance,
#' # then we substitute the values less than 1 for 1
#' HFP_Mexico[HFP_Mexico < 1] <- 1
#' plot(HFP_Mexico)
#'
#' test <- MK_ProtConn(nodes = Protected_areas, region = region,
#'                     area_unit = "ha",
#'                     distance = list(type= "least-cost", resistance = HFP_Mexico),
#'                     distance_thresholds = c(50000, 10000),
#'                     probability = 0.5, transboundary = 50000,
#'                     LA = NULL, plot = TRUE,
#'                     write = NULL, intern = FALSE)
#' test$d50000
#' test$d10000
#' #Using java script
#' #Transform resistance to integer
#' HFP_Mexico <- round(HFP_Mexico)
#'
#' test2 <- MK_ProtConn(nodes = Protected_areas, region = region,
#'                     area_unit = "ha",
#'                     distance = list(type= "least-cost", resistance = HFP_Mexico,
#'                     least_cost.java = TRUE,
#'                     cores.java = 4, ram.java = NULL),
#'                     distance_thresholds = c(50000, 10000),
#'                     probability = 0.5, transboundary = 50000,
#'                     LA = NULL, plot = TRUE,
#'                     write = NULL, intern = FALSE)
#' test2$d50000
#' test2$d10000
#'
#' }
#' @importFrom sf st_as_sf st_cast st_zm st_simplify st_buffer write_sf st_intersection st_difference st_area
#' @importFrom magrittr %>%
#' @importFrom rmapshaper ms_dissolve ms_simplify
#' @importFrom dplyr progress_estimated
#' @importFrom raster raster crop
#' @importFrom fasterize fasterize
#' @importFrom tibble as_tibble
#' @importFrom purrr compact map
#' @importFrom future plan multiprocess availableCores
#' @importFrom furrr future_map
#' @importFrom formattable formattable formatter style color_tile as.htmlwidget
#' @importFrom methods as
#' @importFrom rlang .data
MK_ProtConn <- function(nodes,
                        region,
                        thintersect = NULL,
                        area_unit = "m2",
                        distance = list(type= "centroid", resistance = NULL),
                        distance_thresholds,
                        probability,
                        transboundary = NULL,
                        transboundary_type = "nodes",
                        protconn_bound = FALSE,
                        LA = NULL,
                        geom_simplify = FALSE,
                        plot = FALSE,
                        write = NULL,
                        parallel = NULL,
                        intern = TRUE){
  options(warn = -1)
  err <- tryCatch(detach("package:plyr", unload = TRUE), error = function(err)err)

  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (missing(region)) {
    stop("error missing file of region")
  } else {
    if (is.numeric(region) | is.character(region)) {
      stop("error missing file of region")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(isTRUE(intern)){
    message("Step 1. Reviewing parameters")
  }

  base_param1 <- input_grid(node = nodes, landscape = region, unit = area_unit,
                            bdist = if(is.null(transboundary)){0} else{transboundary},
                            xsimplify = geom_simplify)

  if(class(base_param1)[1] != "input_grid"){
    stop("error in nodes or region shapefile")
  }

  base_param2 <- metric_class(metric = "ProtConn",
                              thintersect = thintersect,
                              distance_threshold = distance_thresholds,
                              probability = probability,
                              transboundary = transboundary,
                              distance = distance)

  if(class(base_param2)[1] != "MK_Metric"){
    stop("error in metric parameters")
  }

  nodes.1 <- tryCatch(Protconn_nodes(x = base_param1@region,
                                     y = base_param1@nodes,
                                     buff = max(transboundary),
                                     method = transboundary_type,
                                     xsimplify = TRUE,
                                     metrunit = base_param1@area_unit,
                                     protconn_bound = protconn_bound), error = function(err)err)
  if(inherits(nodes.1, "error")){
    stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
  }

  if(is.list(nodes.1)){
    if(base_param2@distance$type %in% c("least-cost", "commute-time")){
      if(is.null(base_param2@distance$resistance)){
        stop("error, you need a resistance raster")
      } else {
        centroid <- st_centroid(nodes.1[[1]])
        mask <- st_convex_hull(st_union(centroid)) %>%
          st_buffer(res(base_param2@distance$resistance)[1]*30)
        resist <- crop(base_param2@distance$resistance, as(mask, 'Spatial'))
      }
    } else {
      resist <- NULL
    }
  } else {
    resist <- NULL
  }

  if(is.null(LA)){
    LA <- as.numeric(st_area(base_param1@region)) %>%
      unit_convert(., "m2", base_param1@area_unit)
  }


  base_param3 <- list(base_param1, base_param2, nodes.1, resist, LA)

  if (isTRUE(intern)){
    if(length(base_param3[[2]]@transboundary)>1){
      pb <- dplyr::progress_estimated(length(base_param3[[2]]@transboundary), 0)
      message("Step 2. Processing ProtConn metric. Progress estimated:")
    } else {
      message("Step 2. Processing ProtConn metric")
    }
  }

  if (is.null(parallel)) {
    ProtConn_res <- tryCatch(lapply(base_param3[[2]]@transboundary, function(n){
      if (isTRUE(intern) & length(base_param3[[2]]@transboundary) > 1) {
        pb$tick()$print()
      }

      if(is.list(nodes.1)){
        if(n != max(transboundary)){
          nodes.2 <- tryCatch(Protconn_nodes(x = base_param3[[1]]@region,
                                             y = nodes.1[[1]][which(nodes.1[[1]]$type != "Region"),],
                                             buff = n,
                                             method = transboundary_type,
                                             xsimplify = TRUE,
                                             metrunit = base_param3[[1]]@area_unit,
                                             protconn_bound = protconn_bound), error = function(err)err)
          if(inherits(nodes.2, "error")){
            stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
          }
        } else {
          nodes.2 <- nodes.1
        }

        distance.1 <- tryCatch(protconn_dist(x = nodes.2[[1]], id = "OBJECTID",
                                             y = base_param3[[2]]@distance,
                                             r = base_param3[[1]]@region,
                                             resistance = base_param3[[4]]),
                               error = function(err)err)
        if(inherits(distance.1, "error")){
          stop("error distance. Check topology errors or resistance raster")
        }

          result <- lapply(base_param3[[2]]@distance_threshold, function(d){
            DataProtconn <- get_protconn_grid(x = nodes.2,
                                              y = distance.1,
                                              p = base_param3[[2]]@probability,
                                              pmedian = TRUE,
                                              d = d,
                                              LA = base_param3[[5]], bound = protconn_bound)
            DataProtconn <- round(DataProtconn, 4)

            if(length(which(DataProtconn[5:ncol(DataProtconn)] > 100)) > 0){
              DataProtconn[1, which(DataProtconn[5:ncol(DataProtconn)] > 100) + 5] <- 100
            }

            ##
            DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
            DataProtconn_2$Indicator <- row.names(DataProtconn_2)
            DataProtconn_2$Indicator[1] <- "EC(PC)"
            DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
            DataProtconn_2$Indicator[4] <- "Protected surface"
            DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

            Value <- DataProtconn_2[c(1:4), 1]

            if(Value[1]%%1 == 0){
              Value <- c(formatC(as.numeric(Value[1]), format="d"),
                         formatC(as.numeric(Value[2]), format="e"),
                         formatC(as.numeric(Value[3]), format="d"),
                         formatC(as.numeric(Value[4]), format="d"))
            } else {
              Value <- c(formatC(as.numeric(Value[1]), format="f", digits = 2),
                         formatC(as.numeric(Value[2]), format="e"),
                         formatC(as.numeric(Value[3]), format="f", digits = 2),
                         formatC(as.numeric(Value[4]), format="f", digits = 2))
            }

            DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
            rownames(DataProtconn_2) <- NULL
            #
            DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
            names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
            DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
            rownames(DataProtconn_3) <- NULL

            DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                          list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                               `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                               `Percentage` = color_tile("white", "orange")))
            if(isTRUE(plot)){
              DataProtconn_plot <- plotprotconn(DataProtconn, d)
              result_lista <- list( "Protected Connected (Viewer Panel)" = DataProtconn_4,
                                    "ProtConn Plot" = DataProtconn_plot)
            } else {
              result_lista <- DataProtconn_4
            }
            return(result_lista)
          })

        names(result) <- paste0("d", distance_thresholds)

      } else if(is.numeric(nodes.1)){ #Just exist only one node in the region
        if(n != max(transboundary)){
          nodes.2 <- tryCatch(Protconn_nodes(x = base_param3[[1]]@region,
                                             y = nodes.1[[1]][which(nodes.1[[1]]$type != "Region"),],
                                             buff = n,
                                             xsimplify = TRUE,
                                             metrunit = base_param3[[1]]@area_unit,
                                             protconn_bound = protconn_bound), error = function(err)err)
          if(inherits(nodes.2, "error")){
            stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
          }
        } else {
          nodes.2 <- nodes.1
        }

        result <- lapply(base_param3[[2]]@distance_threshold, function(d){
          DataProtconn <- data.frame(ECA = if(nodes.2 > LA){LA}else{nodes.2},
                                     PC = NA,
                                     LA = LA,
                                     Protected.surface = nodes.2,
                                     Prot = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2/LA)},
                                     Unprotected = if((100 - (100 * (nodes.2 / LA))) < 0){0}else{100 - (100 * (nodes.2 / LA))},
                                     ProtConn = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2 / LA)},
                                     ProtUnconn = 0,
                                     ProtUnconn_Design = 0,
                                     ProtConn_Bound = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2 / LA)},
                                     RelConn = NA,
                                     ProtConn_Prot = NA,
                                     ProtConn_Trans = NA,
                                     ProtConn_Unprot = NA,
                                     ProtConn_Within = NA,
                                     ProtConn_Contig = NA,
                                     ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                     ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)

          DataProtconn[,c(1,3,4:7,10,12)] <- round(DataProtconn[,c(1,3,4:7,10,12)], 4)

          ##
          DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
          DataProtconn_2$Indicator <- row.names(DataProtconn_2)
          DataProtconn_2$Indicator[1] <- "EC(PC)"
          DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
          DataProtconn_2$Indicator[4] <- "Protected surface"
          DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

          Value <- DataProtconn_2[c(1:4),1]

          if(Value[1]%%1 == 0){
            Value <- c(formatC(as.numeric(Value[1]), format="d"),
                       formatC(as.numeric(Value[2]), format="e"),
                       formatC(as.numeric(Value[3]), format="d"),
                       formatC(as.numeric(Value[4]), format="d"))
          } else {
            Value <- c(formatC(as.numeric(Value[1]), format="f", digits = 2),
                       formatC(as.numeric(Value[2]), format="e"),
                       formatC(as.numeric(Value[3]), format="f", digits = 2),
                       formatC(as.numeric(Value[4]), format="f", digits = 2))
          }

          DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
          rownames(DataProtconn_2) <- NULL
          #
          DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
          names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
          DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
          rownames(DataProtconn_3) <- NULL

          DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                        list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `Percentage` = color_tile("white", "orange")))
          return(DataProtconn_4)
        })

        names(result) <- paste0("d", distance_thresholds)
      } else {
        result <- lapply(base_param3[[2]]@distance_threshold, function(d){
          DataProtconn <- data.frame(ECA = NA,
                                     PC = NA,
                                     LA = LA,
                                     Protected.surface = 0,
                                     Prot = 0,
                                     Unprotected = 100,
                                     ProtConn = NA,
                                     ProtUnconn = NA,
                                     ProtUnconn_Design = NA,
                                     ProtConn_Bound = NA,
                                     RelConn = NA,
                                     ProtConn_Prot = NA,
                                     ProtConn_Trans = NA,
                                     ProtConn_Unprot = NA,
                                     ProtConn_Within = NA,
                                     ProtConn_Contig = NA,
                                     ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                     ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)
          ##
          DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
          DataProtconn_2$Indicator <- row.names(DataProtconn_2)
          DataProtconn_2$Indicator[1] <- "EC(PC)"
          DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
          DataProtconn_2$Indicator[4] <- "Protected surface"
          DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

          Value <- DataProtconn_2[c(1:4),1]

          DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
          DataProtconn_2$Value <- round(DataProtconn_2$Value, 4)
          rownames(DataProtconn_2) <- NULL
          #
          DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
          names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
          DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
          rownames(DataProtconn_3) <- NULL

          DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                        list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `Percentage` = color_tile("white", "orange")))
          return(DataProtconn_4)
        })
        names(result) <- paste0("d", distance_thresholds)
      }

      result <- compact(result)

      if(!is.null(write)){
        for (i in 1:length(base_param3[[2]]@distance_threshold)) {
          write.csv(result[[i]][[1]], paste0(write, "_d", base_param3[[2]]@distance_threshold[i],
                                             "_TableProtConn.csv"), row.names = FALSE)
          if(isTRUE(plot)){
            if(!is.character(result[[i]][[2]])){
              tiff(paste0(write, "_d", base_param3[[2]]@distance_threshold[i], '_ProtConn_plot.tif'), width = 806, height = 641)
              print(result[[i]][[2]])
              dev.off()
            }
          }
        }
      }

      return(result)

    }), error = function(err)err)
  } else {
    works <- as.numeric(availableCores())-1
    works <-  if(parallel > works){works}else{parallel}
    plan(strategy = multiprocess, gc = TRUE, workers = works)
    ProtConn_res <- tryCatch(future_map(base_param3[[2]]@transboundary, function(n){
      if(is.list(nodes.1)){
        if(n != max(transboundary)){
          nodes.2 <- tryCatch(Protconn_nodes(x = base_param3[[1]]@region,
                                             y = nodes.1[[1]][which(nodes.1[[1]]$type != "Region"),],
                                             buff = n,
                                             method = transboundary_type,
                                             xsimplify = TRUE,
                                             metrunit = base_param3[[1]]@area_unit,
                                             protconn_bound = protconn_bound), error = function(err)err)
          if(inherits(nodes.2, "error")){
            stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
          }
        } else {
          nodes.2 <- nodes.1
        }

        distance.1 <- tryCatch(protconn_dist(x = nodes.2[[1]], id = "OBJECTID",
                                             y = base_param3[[2]]@distance,
                                             r = base_param3[[1]]@region,
                                             resistance = base_param3[[4]]),
                               error = function(err)err)
        if(inherits(distance.1, "error")){
          stop("error distance. Check topology errors or resistance raster")
        }

        result <- future_map(base_param3[[2]]@distance_threshold, function(d){
          DataProtconn <- get_protconn_grid(x = nodes.2,
                                            y = distance.1,
                                            p = base_param3[[2]]@probability,
                                            pmedian = TRUE,
                                            d = d,
                                            LA = base_param3[[5]], bound = protconn_bound)
          DataProtconn <- round(DataProtconn, 4)

          if(length(which(DataProtconn[5:ncol(DataProtconn)] > 100)) > 0){
            DataProtconn[1,which(DataProtconn[5:ncol(DataProtconn)] > 100) + 5] <- 100
          }
          ##
          DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
          DataProtconn_2$Indicator <- row.names(DataProtconn_2)
          DataProtconn_2$Indicator[1] <- "EC(PC)"
          DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
          DataProtconn_2$Indicator[4] <- "Protected surface"
          DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

          Value <- DataProtconn_2[c(1:4),1]

          if(Value[1]%%1 == 0){
            Value <- c(formatC(as.numeric(Value[1]), format="d"),
                       formatC(as.numeric(Value[2]), format="e"),
                       formatC(as.numeric(Value[3]), format="d"),
                       formatC(as.numeric(Value[4]), format="d"))
          } else {
            Value <- c(formatC(as.numeric(Value[1]), format="f", digits = 2),
                       formatC(as.numeric(Value[2]), format="e"),
                       formatC(as.numeric(Value[3]), format="f", digits = 2),
                       formatC(as.numeric(Value[4]), format="f", digits = 2))
          }

          DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
          rownames(DataProtconn_2) <- NULL
          #
          DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
          names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
          DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
          rownames(DataProtconn_3) <- NULL

          DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                        list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `Percentage` = color_tile("white", "orange")))
          if(isTRUE(plot)){
            DataProtconn_plot <- plotprotconn(DataProtconn, d)
            result_lista <- list( "Protected Connected (Viewer Panel)" = DataProtconn_4,
                                  "ProtConn Plot" = DataProtconn_plot)
          } else {
            result_lista <- DataProtconn_4
          }
          return(result_lista)
        })
        names(result) <- paste0("d", distance_thresholds)

      } else if(is.numeric(nodes.1)){
        if(n != max(transboundary)){
          nodes.2 <- tryCatch(Protconn_nodes(x = base_param3[[1]]@region,
                                             y = nodes.1[[1]][which(nodes.1[[1]]$type != "Region"),],
                                             buff = n,
                                             xsimplify = TRUE,
                                             metrunit = base_param3[[1]]@area_unit,
                                             protconn_bound = protconn_bound), error = function(err)err)
          if(inherits(nodes.2, "error")){
            stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
          }
        } else {
          nodes.2 <- nodes.1
        }
        result <- future_map(base_param3[[2]]@distance_threshold, function(d){
          DataProtconn <- data.frame(ECA = if(nodes.2 > LA){LA}else{nodes.2},
                                     PC = NA,
                                     LA = LA,
                                     Protected.surface = nodes.2,
                                     Prot = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2/LA)},
                                     Unprotected = if((100 - (100 * (nodes.2 / LA))) < 0){0}else{100 - (100 * (nodes.2 / LA))},
                                     ProtConn = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2 / LA)},
                                     ProtUnconn = NA,
                                     ProtUnconn_Design = NA,
                                     ProtConn_Bound = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2 / LA)},
                                     RelConn = NA,
                                     ProtConn_Prot = NA,
                                     ProtConn_Trans = NA,
                                     ProtConn_Unprot = NA,
                                     ProtConn_Within = NA,
                                     ProtConn_Contig = NA,
                                     ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                     ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)

          DataProtconn[,c(1,3,4:7,12)] <- round(DataProtconn[,c(1,3,4:7,12)], 4)

          ##
          DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
          DataProtconn_2$Indicator <- row.names(DataProtconn_2)
          DataProtconn_2$Indicator[1] <- "EC(PC)"
          DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
          DataProtconn_2$Indicator[4] <- "Protected surface"
          DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

          Value <- DataProtconn_2[c(1:4),1]

          if(Value[1]%%1 == 0){
            Value <- c(formatC(as.numeric(Value[1]), format="d"),
                       formatC(as.numeric(Value[2]), format="e"),
                       formatC(as.numeric(Value[3]), format="d"),
                       formatC(as.numeric(Value[4]), format="d"))
          } else {
            Value <- c(formatC(as.numeric(Value[1]), format="f", digits = 2),
                       formatC(as.numeric(Value[2]), format="e"),
                       formatC(as.numeric(Value[3]), format="f", digits = 2),
                       formatC(as.numeric(Value[4]), format="f", digits = 2))
          }

          DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
          rownames(DataProtconn_2) <- NULL
          #
          DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
          names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
          DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
          rownames(DataProtconn_3) <- NULL

          DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                        list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `Percentage` = color_tile("white", "orange")))
          return(DataProtconn_4)
        })

        names(result) <- paste0("d", distance_thresholds)
      } else {
        result <- future_map(base_param3[[2]]@distance_threshold, function(d){
          DataProtconn <- data.frame(ECA = NA,
                                     PC = NA,
                                     LA = LA,
                                     Protected.surface = 0,
                                     Prot = 0,
                                     Unprotected = 100,
                                     ProtConn = NA,
                                     ProtUnconn = NA,
                                     ProtUnconn_Design = NA,
                                     ProtConn_Bound = NA,
                                     RelConn = NA,
                                     ProtConn_Prot = NA,
                                     ProtConn_Trans = NA,
                                     ProtConn_Unprot = NA,
                                     ProtConn_Within = NA,
                                     ProtConn_Contig = NA,
                                     ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                     ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)
          ##
          DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
          DataProtconn_2$Indicator <- row.names(DataProtconn_2)
          DataProtconn_2$Indicator[1] <- "EC(PC)"
          DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
          DataProtconn_2$Indicator[4] <- "Protected surface"
          DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

          Value <- DataProtconn_2[c(1:4),1]

          DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
          DataProtconn_2$Value <- round(DataProtconn_2$Value, 4)
          rownames(DataProtconn_2) <- NULL
          #
          DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
          names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
          DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
          rownames(DataProtconn_3) <- NULL

          DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                        list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                             `Percentage` = color_tile("white", "orange")))
          return(DataProtconn_4)
        })
        names(result) <- paste0("d", distance_thresholds)
      }

      result <- compact(result)

      if(!is.null(write)){
        for (i in 1:length(base_param3[[2]]@distance_threshold)) {
          write.csv(result[[i]][[1]], paste0(write, "_d", base_param3[[2]]@distance_threshold[i],
                                             "_TableProtConn.csv"), row.names = FALSE)
          if(isTRUE(plot)){
            if(!is.character(result[[i]][[2]])){
              tiff(paste0(write, "_d", base_param3[[2]]@distance_threshold[i], '_ProtConn_plot.tif'), width = 806, height = 641)
              print(result[[i]][[2]])
              dev.off()
            }
          }
        }
      }

      return(result)

    }, .progress = intern), error = function(err)err)
    close_multiprocess(works)
    }

  if (inherits(ProtConn_res, "error")){
    stop(ProtConn_res)
  } else {
    if(length(transboundary) > 1){
      names(ProtConn_res) <- paste0("Transboundary_", transboundary)
    } else {
      ProtConn_res <- ProtConn_res[[1]]
      if(length(distance_thresholds) == 1){
          ProtConn_res <- ProtConn_res[[1]]
        }
      }
  }

  return(ProtConn_res)
}
