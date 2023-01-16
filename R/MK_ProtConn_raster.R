#' Protected Connected (ProtConn) using raster nodes
#'
#' Estimate Protected Connected (ProtConn) indicator and fractions for one region.
#' @param nodes raster. The file must have a projected coordinate system.
#' @param region object of class raster, sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
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
#' @param delta logical. Estimate the contribution of each node to the ProtConn value in the region.
#' @param plot logical. Plot the main ProtConn indicators and fractions, default = FALSE.
#' @param write character. Output folder including the output file name without extension, e.g., "C:/ProtConn/Protfiles".
#' @param parallel numeric. Specify the number of cores to use for parallel processing, default = NULL. Parallelize the function using furrr package and multiprocess
#' plan when there are more than ONE transboundary.
#' @param aggregate_raster numeric.Use this parameter if you work with small high-resolution raster (e.g., less than 150 m) or with very large extents. The value of this parameter is used to resample the raster region, perform buffering and speed up the process of transboundary node selection.
#' @param intern logical. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
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
#' data("Protected_areas_raster", package = "Makurhini")
#' data("regions", package = "Makurhini")
#' region <- regions[2,]
#' test <- MK_ProtConn_raster(nodes = Protected_areas_raster, region = region,
#'                            area_unit = "ha",
#'                            distance = list(type= "centroid"),
#'                            distance_thresholds = c(50000, 10000),
#'                            probability = 0.5, transboundary = 50000,
#'                            LA = NULL, plot = FALSE, parallel = NULL,
#'                            protconn_bound=TRUE,
#'                            write = NULL, intern = FALSE)
#' test
#'
#' }
#' @importFrom sf st_buffer write_sf st_area
#' @importFrom magrittr %>%
#' @importFrom raster raster crop buffer clump res res<-
#' @importFrom purrr compact map
#' @importFrom future plan multiprocess availableCores
#' @importFrom furrr future_map
#' @importFrom formattable formattable formatter style color_tile as.htmlwidget
#' @importFrom methods as
#' @importFrom progressr handlers handler_pbcol progressor
#' @importFrom crayon bgWhite white bgCyan
#' @importFrom terra freq rast vect rasterize as.polygons
MK_ProtConn_raster <- function(nodes,
                               region,
                               area_unit = "m2",
                               distance = list(type= "centroid", resistance = NULL),
                               distance_thresholds,
                               probability,
                               transboundary = NULL,
                               transboundary_type = "nodes",
                               protconn_bound = FALSE,
                               LA = NULL,
                               delta = FALSE,
                               plot = FALSE,
                               write = NULL,
                               parallel = NULL,
                               aggregate_raster = NULL,
                               intern = TRUE){
  options(warn = -1)
  . = NULL
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

  if(isTRUE(intern)){
    message("Step 1. Reviewing parameters")
  }

  if(class(region)[1] != "RasterLayer"){
    region.1 <- st_as_sf(region)

    if(isTRUE(protconn_bound)){
      region.1 <- ms_explode(region.1); region.1$id <-  1:nrow(region.1)
    } else {
      region.1$id <-  1
    }

    region.1 <- terra::rasterize(vect(region.1), rast(nodes), "id") %>% raster()

    region.2 <- region.1
    region.2[!is.na(region.2)] <- 1
  } else {
    region.2 <- region
    region.1 <- clump(region.2)
    region.2[!is.na(region.2)] <- 1
  }

  Non_Transboundary_PA <- region.2 * nodes

  v <- raster::values(Non_Transboundary_PA) %>% unique(.)
  v <- v[which(!is.na(v))]

  #if length(v) < 1 = no patches
  if(length(v) < 1){
    loop <- 1:length(distance_thresholds)

    if(is.null(LA)){
      LA <- terra::freq(rast(region.1))[,2:3]

      if(class(LA)[1] == "numeric"){
        LA <- LA[[2]]
      } else {
        LA <- LA[,2]
      }

      LA <- sum(unit_convert(LA * res(nodes)[1]^2,
                             "m2", area_unit))
    }

    if (isTRUE(intern) & length(transboundary) == 1 &
        length(distance_thresholds) > 1) {
      p <- progressor(along = 1:length(distance_thresholds))
    }

    ProtConn_result <- lapply(loop, function(d){
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
      if (isTRUE(intern) & length(transboundary) == 1 &
          length(distance_thresholds) > 1) {
        p()
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
    ProtConn_result <- ProtConn_result[[1]]
    message(paste0("Warning message: No nodes found in the region"))
  }
  else {

    raster_referencia <- if(transboundary_type == "nodes"){Non_Transboundary_PA} else {region.2}

    if(!is.null(aggregate_raster)){
      region.3 <- raster_referencia * 0
      res(region.3) <- aggregate_raster
      region.3 <- resample(raster_referencia, region.3, method = "ngb") %>%
        raster::buffer(., transboundary)

      r <- raster(extent(region.3), res = res(raster_referencia))
      region.3 <- resample(region.3, r, method = "ngb")
      crs(region.3) <- crs(raster_referencia)

    } else {
      region.3 <- as.polygons(rast(raster_referencia)) %>%
        st_as_sf() %>% ms_dissolve() %>% st_buffer(., 0)

      region.3 <- st_buffer(region.3, transboundary)

      region.3 <- terra::rasterize(vect(region.3), rast(raster_referencia)) %>% raster(.)
    }

    region.4 <- reclassify(region.2, c(NA,-Inf,1, -Inf,1,NA))
    region.4 <- region.3 * region.4

    Transboundary_PA <- region.4 * nodes

    ##Table of nodes and attributes
    Val <- terra::freq(rast(Non_Transboundary_PA))[,2:3]

    if(class(Val)[1] == "numeric"){
      Val_id <- Val[[1]]
      Val_nodes <- Val[[2]]
    } else {
      Val_id <- Val[,1]
      Val_nodes <- Val[,2]
    }

    if(length(Val_nodes) > 1){
      nodes.1 <- data.frame(id = 1:length(Val_nodes),
                            id2 = Val_id,
                            attribute = unit_convert(Val_nodes * res(nodes)[1]^2,
                                                     "m2", area_unit),
                            type = "Non-Transboundary")
      #centroid
      nodes_centroid <- rast_clump_points(x = Non_Transboundary_PA,
                                          parallel = parallel,
                                          centroid_geometry = TRUE)
      nodes_centroid$id <- 1
      nodes_centroid$id2 <- nodes.1$id2
      nodes_centroid$type <- "Non-Transboundary"

      ##
      Val <- terra::freq(rast(Transboundary_PA))[,2:3]

      if(class(Val)[1] == "numeric"){
        Val_id <- Val[[1]]
        Val_nodes <- Val[[2]]
      } else {
        Val_id <- Val[,1]
        Val_nodes <- Val[,2]
      }

      if(length(Val_nodes) > 1){
        Val_nodes <- data.frame(id = 1:length(Val_nodes),
                                id2 = Val_id,
                                attribute = unit_convert(Val_nodes * res(nodes)[1]^2,
                                                         "m2", area_unit),
                                type = "Transboundary")
        Val_nodes$id <- Val_nodes$id + max(nodes.1$id)
        nodes.1 <- rbind(nodes.1, Val_nodes)
        #
        nodes_centroid.2 <- rast_clump_points(x = Transboundary_PA,
                                              parallel = parallel,
                                              centroid_geometry = TRUE)
        nodes_centroid.2$id <- 1
        nodes_centroid.2$id2 <- Val_nodes[,1]
        nodes_centroid.2$type <- "Transboundary"
        nodes_centroid <- rbind(nodes_centroid, nodes_centroid.2)
        nodes_centroid$id <- 1:nrow(nodes_centroid)
      }

      if(isTRUE(protconn_bound)){
        Val <- terra::freq(rast(region.1))[,2:3]

        if(class(Val)[1] == "numeric"){
          Val_id <- Val[[1]]
          Val_nodes <- Val[[2]]
        } else {
          Val_id <- Val[,1]
          Val_nodes <- Val[,2]
        }

        Val_nodes <- data.frame(id = 1:length(Val_nodes),
                                id2 = NA,
                                attribute = unit_convert(Val_nodes * res(nodes)[1]^2,
                                                         "m2", area_unit),
                                type = "Region")
        Val_nodes$id <- Val_nodes$id + max(nodes.1$id)
        nodes.1 <- rbind(nodes.1, Val_nodes)

        region_centroid <- rast_clump_points(x = region.1,
                                             parallel = parallel,
                                             centroid_geometry = TRUE)
        region_centroid$id <- 1
        region_centroid$id2 <- NA
        region_centroid$type <- "Region"

        nodes_centroid <- rbind(nodes_centroid, region_centroid)
        nodes_centroid$id <- 1:nrow(nodes_centroid)
      }

    } else {
      if(length(Val_nodes) == 1){
        nodes.1 <- data.frame(id = 1,
                              id2 = Val_id,
                              attribute = unit_convert(Val_nodes * res(nodes)[1]^2,
                                                       "m2", area_unit),
                              type = "Non-Transboundary")
        nodes_centroid <- NULL
      } else {
        nodes.1 <- "number of nodes less than 2"
        nodes_centroid <- "number of nodes less than 2"
      }
    }

    if(!is.character(nodes.1)){
      if(nrow(nodes.1) > 1){
        if(distance$type %in% c("least-cost", "commute-time")){
          if(is.null(distance$resistance)){
            stop("error, you need a resistance raster")
          } else {
            mask <- st_convex_hull(st_union(nodes_centroid)) %>%
              st_buffer(res(distance$resistance)[1]*30)
            resist <- crop(distance$resistance, as(mask, 'Spatial'))
          }
        } else {
          resist <- NULL
        }
      } else {
        resist <- NULL
      }
    } else {
      resist <- NULL
    }

    if(is.null(LA)){
      LA <- terra::freq(rast(region.1))[,2:3]

      if(class(LA)[1] == "numeric"){
        LA <- LA[[2]]
      } else {
        LA <- LA[,2]
      }

      LA <- sum(unit_convert(LA * res(nodes)[1]^2,
                             "m2", area_unit))

      if(!is.character(nodes.1) & nrow(nodes.1) == 1){
        if(nodes.1$attribute >= LA){
          nodes.1$attribute <- LA
        }
      }
    }

    if (isTRUE(intern)){
      handlers(global = T, append = TRUE)
      handlers(handler_pbcol(complete = function(s) crayon::bgYellow(crayon::white(s)),
                             incomplete = function(s) crayon::bgWhite(crayon::black(s)),
                             intrusiveness = 2))
      if(length(transboundary) > 1 | length(distance_thresholds) > 1){
        message("Step 2. Processing ProtConn metric. Progress estimated:")
      } else {
        message("Step 2. Processing ProtConn metric")
      }
    }

    ProtConn_Estimation <- function(patches = NULL, n = NULL,
                                    distance = NULL, intern = TRUE,
                                    distance_thresholds = NULL,
                                    probability = NULL,
                                    delta = NULL,
                                    LA = NULL, works = NULL,
                                    write = NULL){
      if(class(patches[[1]])[1] == "sf"){
        patches[[1]]$attribute <- patches[[2]]$attribute
        nodes.2 <- patches[[1]]

        distance.1 <- tryCatch(protconn_dist(x = nodes.2, id = "id",
                                             y = distance,
                                             r = region,
                                             resistance = resist),
                               error = function(err)err)

        if(inherits(distance.1, "error")){
          stop("error distance. Check topology errors or resistance raster")
        }

        if (isTRUE(intern) & length(transboundary) == 1 &
            length(distance_thresholds) > 1) {
          p <- progressor(along = 1:length(distance_thresholds))
        }

        loop <- 1:length(distance_thresholds)
        ##
        result <- lapply(loop, function(d){
          d.2 <- distance_thresholds[d]
          DataProtconn <- get_protconn_grid(x = list(patches[[1]], patches[[1]]$attribute),
                                            y = distance.1,
                                            p = probability,
                                            pmedian = TRUE,
                                            d = d.2,
                                            LA = LA, bound = protconn_bound)
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

          if (isTRUE(intern) & length(transboundary) == 1 &
              length(distance_thresholds) > 1) {
            p()
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
            DataProtconn_plot <- plotprotconn(DataProtconn, d.2)
            result_lista <- list( "Protected Connected (Viewer Panel)" = DataProtconn_4,
                                  "ProtConn Plot" = DataProtconn_plot)
          } else {
            result_lista <- DataProtconn_4
          }
          return(result_lista)
        })

      }
      else if (is.null(patches[[1]])) {
        nodes.2 <- patches[[2]]$attribute

        if(isTRUE(intern) & length(transboundary) == 1 &
           length(distance_thresholds) > 1) {
          p <- progressor(along = 1:length(distance_thresholds))
        }

        loop <- 1:length(distance_thresholds)

        result <- lapply(loop, function(d){
          d.2 <- distance_thresholds[d]
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
          if (isTRUE(intern) & length(transboundary) == 1 &
              length(distance_thresholds) > 1) {
            p()
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

      }
      else {
        nodes.2 <- patches[[2]]
        if (isTRUE(intern) & length(transboundary) == 1 &
            length(distance_thresholds) > 1) {
          p <- progressor(along = 1:length(distance_thresholds))
        }

        loop <- 1:length(distance_thresholds)

        result <- lapply(loop, function(d){
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
          if (isTRUE(intern) & length(transboundary) == 1 &
              length(distance_thresholds) > 1) {
            p()
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
        names(result) <- paste0("d", distance_thresholds)
        message(paste0("Warning message: No nodes found in the region, transboundary "), n)
      }

      if(isTRUE(delta)){
        if (isTRUE(intern)){
          message("Step 3. Processing Delta ProtConn")
        } else {
          message("Processing Delta ProtConn")
        }

        if(is.character(nodes.2)){
          message("Analysis cannot be completed, no nodes in the region, transboundary ", n)
          deltaProtConn <- paste0("Analysis cannot be completed, no nodes in the region, transboundary ", n)
        } else {
          deltaProtConn <- tryCatch(delta_ProtConn_raster(x= patches[[1]],
                                                          y= patches[[3]],
                                                          distance_thresholds = distance_thresholds,
                                                          probability = probability,
                                                          LA = LA,
                                                          distance = distance,
                                                          resist = resist,
                                                          works = works),
                                    error = function(err)err)
          if(inherits(deltaProtConn, "error")){
            message("Analysis cannot be completed, no nodes in the region, transboundary ", n)
            deltaProtConn <- paste0("Analysis cannot be completed, no nodes in the region, transboundary ", n)
          }
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

      names(result) <- paste0("d", distance_thresholds)
      result <- purrr::compact(result)

      if(!is.null(write)){
        for (i in 1:length(distance_thresholds)) {
          write.csv(result[[i]][[1]], paste0(write, "_d", distance_thresholds[i],
                                             "_TableProtConn.csv"), row.names = FALSE)
          if(isTRUE(plot)){
            if(!is.character(result[[i]][[2]])){
              tiff(paste0(write, "_d", distance_thresholds[i], '_ProtConn_plot.tif'), width = 806, height = 641)
              print(result[[i]][[2]])
              dev.off()
            }
          }
        }
      }
      return(result)
    }

    if(!is.null(parallel)){
      works <- as.numeric(availableCores())-1
      works <-  if(parallel > works){works}else{parallel}
    } else {
      works <- NULL
    }

    ProtConn_Estimation_progress <- function(xs) {
      if (isTRUE(intern) & length(transboundary) > 1) {
        p <- progressor(along = xs)
      }

      y <- lapply(xs, function(x){
        x.1 <- ProtConn_Estimation(patches = list(nodes_centroid, nodes.1, nodes),
                                   n = transboundary[x],
                                   distance = distance,
                                   distance_thresholds = distance_thresholds,
                                   probability = probability,
                                   delta = delta,
                                   LA = LA, works = works,
                                   write = write, intern = intern)
        if (isTRUE(intern) & length(transboundary) > 1) {
          p()
        }
        return(x.1)
      })
      close_multiprocess(works)
      return(y)
    }

    ProtConn_res <- tryCatch(ProtConn_Estimation_progress(xs = 1:length(transboundary)),
                             error = function(err)err)

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
  }

  return(ProtConn_res)
}
