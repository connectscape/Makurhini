#' Protected Connected (ProtConn)
#'
#' @description
#' This function calculates the Protected Connected indicator (ProtConn) for a region, its fractions and the importance (contribution) of each protected area to maintain connectivity in the region under one or more distance thresholds.
#' @param nodes object of class \code{sf, sfc, sfg, spatialPolygonsDataFrame}. Spatial data of vector type that normally contains the spatial limits of protected areas. It must be in a projected coordinate system.
#' @param region object of class \code{sf, sfc, sfg, spatialPolygonsDataFrame}. Polygon delimiting the region or study area. It must be in a projected coordinate system.
#' @param area_unit \code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.
#' @param distance A \code{list} of parameters to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. \cr
#'  This list must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. For example: \code{distance(type = "least-cost", resistance = raster_resistance)}. \cr
#' To see more arguments see the \link[Makurhini]{distancefile} function.
#' @param distance_thresholds A \code{numeric} indicating the dispersal distance or distances (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range. Can be the same length as the \code{distance_thresholds} parameter.
#' @param probability A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the \code{"PC"} metric. If \code{probability = NULL}, then a probability of 0.5 will be used.
#' @param transboundary \code{numeric}. Buffer to select polygons (e.g., PAs) in a second round. The selected polygons will have an attribute value = 0, i.e., their contribution for connectivity would be as stepping stones (Saura et al. 2017). One cross-border value or one for each threshold distance can be set.
#' @param transboundary_type \code{character}. Two options: \code{"nodes" (methodology from Saura et al. 2017)} or \code{"region"}.\cr
#' - If it is \code{"nodes"}, the transboundary is built from the limits of the nodes present in the region (default).
#' - If it is \code{"region"}, is selected the transboundary is built from the limits of the region.
#' @param protconn_bound \code{logical}. If \code{TRUE} then the fractions ProtUnConn[design] and ProtConn[bound] will be estimated.
#' @param geom_simplify \code{logical}. Slightly simplify the region and nodes geometries.
#' @param delta \code{logical}. Estimate the contribution of each node to the ProtConn value in the region.
#' @param plot \code{logical}. Plot the main ProtConn indicators and fractions, default = \code{FALSE}.
#' @param write \code{character}. Output folder including the output file name without extension, e.g., \code{"C:/ProtConn/Protfiles"}.
#' @param parallel \code{numeric}. Specify the number of cores to use for parallel processing, default = NULL. Parallelize the function using furrr package and multiprocess
#' plan when there are more than one transboundary.
#' @param intern \code{logical}. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @return
#' Table with the following ProtConn values: \code{ECA, Prot, ProtConn, ProtUnconn, RelConn, ProtUnConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within],
#'  ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot], ProtConn[Within][land], ProtConn[Contig][land],
#'  ProtConn[Unprot][land], ProtConn[Trans][land]} \cr
#' \cr
#' - If plot \bold{is not NULL} a \code{list} is returned with the ProtConn table and a plots.
#' - If \bold{delta} is \code{TRUE} then it returns an sf class object with the importance value (contribution to ProtConn) for each node in the region.
#' @references
#' Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.\cr
#' Saura, S., Bertzky, B., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2018). Protected area connectivity: Shortfalls in global targets and country-level priorities. Biological Conservation, 219(October 2017), 53–67.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(sf)
#'
#' data("Protected_areas", package = "Makurhini")
#' data("regions", package = "Makurhini")
#' region <- regions[2,]
#'
#' test <- MK_ProtConn(nodes = Protected_areas, region = region,
#'                     area_unit = "ha",
#'                     distance = list(type= "centroid"),
#'                     distance_thresholds = c(50000, 10000),
#'                     probability = 0.5, transboundary = 50000,
#'                     plot = TRUE, parallel = NULL,
#'                     protconn_bound = TRUE,
#'                     delta = TRUE,
#'                     write = NULL, intern = TRUE)
#' test
#'
#' #Least-cost distances using a human foot print
#' library(raster)
#' HFP_Mexico <- raster(system.file("extdata", "HFP_Mexico.tif",
#'                     package = "Makurhini", mustWork = TRUE))
#' mask_1 <- as(extent(Protected_areas), 'SpatialPolygons')
#' crs(mask_1) <- crs(Protected_areas)
#' mask_1 <- buffer(mask_1, 20000)
#' HFP_Mexico <- crop(HFP_Mexico, mask_1)
#' #If least_cost.java is TRUE, then resistance must bee an integer raster (i.e., integer values).
#' HFP_Mexico <- round(HFP_Mexico)
#'
#' test2 <- MK_ProtConn(nodes = Protected_areas, region = region,
#'                     area_unit = "ha",
#'                     distance = list(type = "least-cost", resistance = HFP_Mexico,
#'                     least_cost.java = TRUE,
#'                     cores.java = 4, ram.java = NULL),
#'                     distance_thresholds = c(50000, 10000),
#'                     probability = 0.5, transboundary = 50000,
#'                     plot = TRUE,
#'                     write = NULL, intern = FALSE)
#' test2$d50000
#' test2$d10000
#'
#' }
#' @importFrom sf st_buffer write_sf st_area
#' @importFrom magrittr %>%
#' @importFrom raster raster crop
#' @importFrom purrr compact map
#' @importFrom future plan multicore multisession availableCores
#' @importFrom furrr future_map
#' @importFrom formattable formattable formatter style color_tile as.htmlwidget
#' @importFrom methods as
#' @importFrom utils installed.packages txtProgressBar setTxtProgressBar
MK_ProtConn <- function(nodes = NULL,
                        region = NULL,
                        area_unit = "m2",
                        distance = list(type= "centroid", resistance = NULL),
                        distance_thresholds = NULL,
                        probability = 0.5,
                        transboundary = NULL,
                        transboundary_type = "nodes",
                        protconn_bound = FALSE,
                        geom_simplify = FALSE,
                        delta = FALSE,
                        plot = FALSE,
                        write = NULL,
                        parallel = NULL,
                        intern = TRUE){
  options(warn = -1); . = NULL
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

  if(isFALSE(parallel)){
    parallel <- NULL
  } else if (isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
  } else if((!is.null(parallel))){
    if(!is.numeric(parallel)){
      stop("if you use parallel argument then you need a numeric value")
    }
  } else {
    parallel <- NULL
  }

  if(isTRUE(intern)){
    message("Step 1. Reviewing parameters")
  }

  base_param1 <- tryCatch(input_grid(node = nodes, landscape = region, unit = area_unit,
                                     bdist = if(is.null(transboundary)){0} else{transboundary},
                                     xsimplify = geom_simplify), error = function(err)err)

  if(inherits(base_param1, "error")){
    stop("error in nodes or region shapefile")
  }

  base_param2 <- metric_class(metric = "ProtConn",
                              distance_threshold = distance_thresholds,
                              probability = probability,
                              transboundary = transboundary,
                              distance = distance)

  if(nrow(base_param1@nodes) > 0){
    nodes.1 <- tryCatch(Protconn_nodes(x = base_param1@region,
                                       y = base_param1@nodes,
                                       buff = max(transboundary),
                                       method = transboundary_type,
                                       xsimplify = geom_simplify,
                                       metrunit = base_param1@area_unit,
                                       protconn_bound = protconn_bound,
                                       delta = delta), error = function(err)err)

    if(inherits(nodes.1, "error")){
      stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
    }

    if(is.numeric(nodes.1)){
      nodes.delta <- over_poly(x = base_param1@nodes, y = base_param1@region, geometry = TRUE)
    }
  } else {
    nodes.1 <- "No nodes"
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

  LA <- as.numeric(st_area(base_param1@region)) %>%
    unit_convert(., "m2", base_param1@area_unit)

  if(is.numeric(nodes.1)){
    if(nodes.1 >= LA){
      nodes.1 <- LA
    }
  }

  base_param3 <- list(base_param1, base_param2, nodes.1, resist, LA)

  if (isTRUE(intern)){
    if(length(base_param3[[2]]@transboundary) >1 | length(base_param3[[2]]@distance_threshold) > 1){
      message("Step 2. Processing ProtConn metric. Progress estimated:")
    } else {
      message("Step 2. Processing ProtConn metric")
    }
  }

  if (isTRUE(intern) & length(base_param3[[2]]@transboundary) > 1) {
    pb <- txtProgressBar(0,length(base_param3[[2]]@transboundary), style = 3)
  }

  if(is.null(parallel)){
    ProtConn_res <- tryCatch(lapply(1:length(base_param3[[2]]@transboundary), function(x){
      x.1 <- ProtConn_Estimation(base_param3 = base_param3,
                                 base_param1 = base_param1,
                                 nodes.delta = nodes.delta,
                                 n = base_param3[[2]]@transboundary[x],
                                 delta = delta,
                                 transboundary_type = transboundary_type,
                                 transboundary = transboundary,
                                 protconn_bound = protconn_bound,
                                 write = write,
                                 LA = LA,
                                 plot = plot,
                                 intern = intern)

      if (isTRUE(intern) & length(base_param3[[2]]@transboundary) > 1) {
        setTxtProgressBar(pb, x)
      }
      return(x.1)
    }), error = function(err) err)
  } else {
    works <- as.numeric(availableCores())-1; works <-  if(parallel > works){works}else{parallel}
    if(.Platform$OS.type == "unix") {
      strat <- future::multicore
    } else {
      strat <- future::multisession
    }
    plan(strategy = strat, gc = TRUE, workers = works)
    ProtConn_res <- tryCatch(future_map(1:length(base_param3[[2]]@transboundary), function(x){
      x.1 <- ProtConn_Estimation(base_param3 = base_param3,
                                 base_param1 = base_param1,
                                 nodes.delta = nodes.delta,
                                 n = base_param3[[2]]@transboundary[x],
                                 delta = delta,
                                 transboundary_type = transboundary_type,
                                 transboundary = transboundary,
                                 protconn_bound = protconn_bound,
                                 LA = LA,
                                 plot = plot,
                                 intern = intern,
                                 write = write)
      return(x.1)
    }, .progress = intern), error = function(err) err)
    close_multiprocess(works)
  }

  if(inherits(ProtConn_res, "error")){
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

  if(isTRUE(intern)){
    message("Done!")
  }

  return(ProtConn_res)
}

#' Protconn estimation
#'
#' @param base_param3 list
#' @param base_param1 list
#' @param n numeric
#' @param delta logical
#' @param nodes.delta sf
#' @param transboundary_type character
#' @param transboundary numeric
#' @param protconn_bound logical
#' @param LA numeric
#' @param plot logical
#' @param write character
#' @param intern logical
#' @importFrom purrr compact
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom grDevices tiff dev.off
#' @importFrom formattable formattable formatter style color_tile as.htmlwidget
#' @keywords internal

ProtConn_Estimation <- function(base_param3 = NULL,
                                base_param1 = NULL,
                                nodes.delta = NULL,
                                n = NULL, delta = FALSE,
                                transboundary_type = NULL,
                                transboundary = NULL,
                                protconn_bound = FALSE,
                                plot = TRUE,
                                LA = NULL,
                                intern = TRUE, write = NULL){
  if(is.list(base_param3[[3]])){
    if(n != max(transboundary)){
      nodes.2 <- tryCatch(Protconn_nodes(x = base_param3[[1]]@region,
                                         y = base_param1@nodes,
                                         buff = n,
                                         method = transboundary_type,
                                         xsimplify = TRUE,
                                         metrunit = base_param3[[1]]@area_unit,
                                         protconn_bound = protconn_bound,
                                         delta = delta), error = function(err)err)
      if(inherits(nodes.2, "error")){
        stop(paste0("error first nodes selection, please check topology errors and you could simplify polygon"))
      }
    } else {
      nodes.2 <- base_param3[[3]]
    }

    distance.1 <- tryCatch(protconn_dist(x = nodes.2[[1]], id = "OBJECTID",
                                         y = base_param3[[2]]@distance,
                                         r = base_param3[[1]]@region,
                                         resistance = base_param3[[4]]),
                           error = function(err)err)

    if(inherits(distance.1, "error")){
      stop("error distance. Check topology errors or resistance raster")
    }

    if(isTRUE(intern) & length(base_param3[[2]]@transboundary) == 1 &
       length(base_param3[[2]]@distance_threshold ) > 1) {
      pb <- txtProgressBar(0,length(base_param3[[2]]@distance_threshold), style = 3)
    }
    result <- lapply(1:length(base_param3[[2]]@distance_threshold), function(x){
      d.2 <- base_param3[[2]]@distance_threshold[x]
      DataProtconn <- get_protconn_grid(x = nodes.2,
                                        y = distance.1,
                                        p = base_param3[[2]]@probability,
                                        pmedian = TRUE,
                                        d = d.2,
                                        LA = base_param3[[5]],
                                        bound = protconn_bound)
      DataProtconn <- round(DataProtconn, 4)

      if(length(which(DataProtconn[5:ncol(DataProtconn)] > 100)) > 0){
        DataProtconn[1, which(DataProtconn[5:ncol(DataProtconn)] > 100) + 5] <- 100
      }

      if(length(which(DataProtconn[5:ncol(DataProtconn)] < 0)) > 0){
        DataProtconn[1, which(DataProtconn[5:ncol(DataProtconn)] < 0) + 5] <- 0
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

      if (isTRUE(intern) & length(base_param3[[2]]@transboundary) == 1 &
          length(base_param3[[2]]@distance_threshold ) > 1) {
        setTxtProgressBar(pb, x)
      }

      DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
      rownames(DataProtconn_2) <- NULL
      DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
      names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
      DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
      rownames(DataProtconn_3) <- NULL

      DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                    list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `Percentage` = color_tile("#FFF3DD", "orange")))

      if(isTRUE(plot)){
        DataProtconn_plot <- tryCatch(plotprotconn(DataProtconn, d = d.2), error = function(err)err)

        if(inherits(DataProtconn_plot, "error")){
          message("The plot could not be performed, check that the ggplot2 and ggpubr packages are installed or updated")
          result_lista <- DataProtconn_4; plot = FALSE
        } else {
          result_lista <- list( "Protected Connected (Viewer Panel)" = DataProtconn_4,
                                "ProtConn Plot" = DataProtconn_plot)
        }
      } else {
        result_lista <- DataProtconn_4
      }
      return(result_lista)
    })
  } else if(is.numeric(base_param3[[3]])) { #Just exist only one node in the region
    nodes.2 <- base_param3[[3]]
    if (isTRUE(intern) & length(base_param3[[2]]@transboundary) == 1 &
        length(base_param3[[2]]@distance_threshold ) > 1) {
      pb <- txtProgressBar(0,length(base_param3[[2]]@distance_threshold), style = 3)
    }

    result <- lapply(1:length(base_param3[[2]]@distance_threshold), function(x){
      d.2 <- base_param3[[2]]@distance_threshold[x]
      DataProtconn <- data.frame(ECA = if(nodes.2 >= LA){LA}else{nodes.2},
                                 PC = if(nodes.2 >= LA){1}else{nodes.2/LA^2},
                                 LA = LA,
                                 Protected.surface = nodes.2,
                                 Prot = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2/LA)},
                                 Unprotected = if((100 - (100 * (nodes.2 / LA))) < 0){0}else{100 - (100 * (nodes.2 / LA))},
                                 ProtConn = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2 / LA)},
                                 ProtUnconn = 0,
                                 ProtUnconn_Design = 0,
                                 ProtConn_Bound = if((100 * (nodes.2 / LA)) > 100){100}else{100 * (nodes.2 / LA)},
                                 RelConn = NA,
                                 ProtConn_Prot = 100,
                                 ProtConn_Trans = NA,
                                 ProtConn_Unprot = NA,
                                 ProtConn_Within = 100,
                                 ProtConn_Contig = NA,
                                 ProtConn_Within_land = NA,
                                 ProtConn_Contig_land = NA,
                                 ProtConn_Unprot_land = NA,
                                 ProtConn_Trans_land = NA)

      DataProtconn[,c(1,3,4:7,10,12)] <- round(DataProtconn[,c(1,3,4:7,10,12)], 4)

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
      if (isTRUE(intern) & length(base_param3[[2]]@transboundary) == 1 &
          length(base_param3[[2]]@distance_threshold ) > 1) {
        setTxtProgressBar(pb, x)
      }
      DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
      rownames(DataProtconn_2) <- NULL
      #
      DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
      names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
      DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
      rownames(DataProtconn_3) <- NULL

      if(isFALSE(protconn_bound)){
        DataProtconn_3 <- DataProtconn_3[-which(DataProtconn_3[,3] %in% c("ProtUnconn_Design", "ProtConn_Bound")),]
      }

      DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                    list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `Percentage` = color_tile("#FFF3DD", "orange")))
      if(DataProtconn_4[[2]][3] == "NA"){
        DataProtconn_4[[2]][3] <- paste(round(LA, 2))
      }
      return(DataProtconn_4)
    })
  } else {
    nodes.2 <- base_param3[[3]]
    if (isTRUE(intern) & length(base_param3[[2]]@transboundary) == 1 &
        length(base_param3[[2]]@distance_threshold ) > 1) {
      pb <- txtProgressBar(0,length(base_param3[[2]]@distance_threshold), style = 3)
    }

    result <- lapply(1:length(base_param3[[2]]@distance_threshold), function(x){
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
      DataProtconn[,3] <- round(DataProtconn[,3], 4)
      DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
      DataProtconn_2$Indicator <- row.names(DataProtconn_2)
      DataProtconn_2$Indicator[1] <- "EC(PC)"
      DataProtconn_2$Indicator[3] <- "Maximum landscape attribute"
      DataProtconn_2$Indicator[4] <- "Protected surface"
      DataProtconn_2$Index <- rep(DataProtconn_2[c(1:4),2], 8)[1:nrow(DataProtconn_2)]

      Value <- DataProtconn_2[c(1:4),1]

      if(Value[3]%%1 == 0){
        Value <- c(formatC(as.numeric(Value[1]), format="d"),
                   formatC(as.numeric(Value[2]), format="e"),
                   formatC(as.numeric(Value[3]), format="d"),
                   formatC(as.numeric(Value[4]), format="d"))
      } else {
        Value <- c(formatC(as.numeric(Value[1]), format="f"),
                   formatC(as.numeric(Value[2]), format="e"),
                   formatC(as.numeric(Value[3]), format="f", digits = 2),
                   formatC(as.numeric(Value[4]), format="f", digits = 0))
      }
      if (isTRUE(intern) & length(base_param3[[2]]@transboundary) == 1 &
          length(base_param3[[2]]@distance_threshold ) > 1) {
        setTxtProgressBar(pb, x)
      }
      DataProtconn_2$Value <- rep(Value, 8)[1:nrow(DataProtconn_2)]
      rownames(DataProtconn_2) <- NULL
      #
      DataProtconn_3 <- DataProtconn_2[5:nrow(DataProtconn_2),c(3:4, 2, 1)]
      names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
      DataProtconn_3[5:nrow(DataProtconn_3), 1:2] <- " "
      rownames(DataProtconn_3) <- NULL

      if(isFALSE(protconn_bound)){
        DataProtconn_3 <- DataProtconn_3[-which(DataProtconn_3[,3] %in% c("ProtUnconn_Design", "ProtConn_Bound")),]
      }

      DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                    list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `Percentage` = color_tile("#FFF3DD", "orange")))

      if(DataProtconn_4[[2]][3] == "NA"){
        DataProtconn_4[[2]][3] <- paste(round(LA,2))
      }

      return(DataProtconn_4)
    })
    names(result) <- paste0("d", base_param3[[2]]@distance_threshold)
    message(paste0("Warning message: No nodes found in the region, transboundary "), n)
  }

  if(isTRUE(delta)){
    if (isTRUE(intern)){
      message("Step 3. Processing Delta ProtConn")
    } else {
      message("Processing Delta ProtConn")
    }

    if(is.character(nodes.2)){
      deltaProtConn <- message(paste0("Analysis cannot be completed, no nodes in the region, transboundary "), n)
    } else {
      deltaProtConn <- delta_ProtConn(x= if(is.numeric(nodes.2)){nodes.delta} else {nodes.2$delta},
                                      y= if(is.numeric(nodes.2)){NULL} else {nodes.2$nodes_diss},
                                      base_param3)
    }

    if(isTRUE(plot) & !is.numeric(nodes.2) & !is.character(nodes.2)){
      for(i in 1:length(result)){
        result[[i]]$'ProtConn_Delta' <- deltaProtConn[[i]]
      }
    } else {
      result <- lapply(1:length(result), function(l){
        l.1 <- list("Protected Connected (Viewer Panel)" = result[[l]],
                    "ProtConn_Delta" = if(is.character(nodes.2)){deltaProtConn} else{deltaProtConn[[l]]})
        return(l.1) })
    }
  }

  names(result) <- paste0("d", base_param3[[2]]@distance_threshold)
  result <- purrr::compact(result)

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
}
