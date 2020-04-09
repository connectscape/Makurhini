#' Multiple Protected Connected (ProtConn)
#'
#' Use the CONEFOR command line to estimate Protected Connected (ProtConn) indexes for multiple ecoregions
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. Protected areas (PAs) shapefile.
#' @param regions object of class sf, sfc, sfg or SpatialPolygons. Ecoregions shapefile.
#' @param thintersect numeric.Threshold of intersection in percentage allowed to select or not a target geometry. Default = 90, if intersection >=90 percentage, the geometry will be selected.
#' @param attribute character. Select the nodes attribute: "Area" = Complete Protected areas; "Intersected area" = Intersected Protected areas; or another specific column name with the nodes attribute, ideally this attribute mus be an area-weighted index, otherwise the interpretation of the protconn index may change.
#' @param area_unit character. Attribute area units. You can set an area unit, "Makurhini::unit_covert()" compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to hectares "ha".
#' @param res_attribute numeric. If the attribute is no equal to "Area" or "Intersected area" then nodes will be converted to raster to extract values in one  process step, you can set the raster resolution, default = 150.
#' @param distance list.See distancefile(). E.g.: list(type = "centroid", resistance = NULL, geometry_out = NULL).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_threshold = c(30000, 50000); sequence distances: distance_threshold = seq(10000,100000, 10000).
#' @param transboundary numeric. Buffer to select polygones in a second round, their attribute value = 0, see Saura et al. 2017.
#' @param keep numeric. Simplification of the region's geometry to speed the select by location ("Makurhini::MK_selectbyloc").
#'  The value can range from 0 to 1 and is the proportion of points to retain (default 0.02). The higher the value,
#'   the higher the speed but the greater uncertainty.
#' @param CI character. A character vector representing the type of confidence intervals that will be estimated. The value should be any subset of the values c("norm","basic", "stud", "perc", "bca") or "all" which will compute all five types of intervals (see, boot::boot.ci())
#' @param plot logical. Plot the main ProtConn indicators and fractions with their standard deviation, default = FALSE.
#' @param write character. Output folder including the output file name without extension, e.g., "C:/ProtConn/Protfiles".
#' @param intern logical. Show the progress of the process, default = TRUE.
#' @param parallel logical. Parallelize the function using furrr package and multiprocess plan, default = FALSE.
#' @return Table with the following ProtConn values: ECA, Prot, ProtConn, ProtUnconn, RelConn, ProtConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within],
#' ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot], ProtConn[Within][land], ProtConn[Contig][land],
#' ProtConn[Unprot][land], ProtConn[Trans][land] \cr
#' \cr
#' *Each indicator is accompanied by six dispersion statistics: standard deviation, standard error and four confidence intervals obtained with a bootstrap-type resampling (Carpenter & Bithell 2000; see, \strong{ProtConnStat()})\cr
#' \cr
#' *If plot is not NULL a list is returned with the ProtConn table and a plots.
#' @references Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S., Bertzky, B., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2018). Protected area connectivity: Shortfalls in global targets and country-level priorities. Biological Conservation, 219(October 2017), 53–67.
#' @export
#' @importFrom magrittr %>%
#' @import sf
#' @import ggplot2
#' @importFrom dplyr progress_estimated
#' @importFrom future multiprocess plan
#' @importFrom furrr future_map
#' @importFrom formattable color_tile area style formattable formatter
#' @importFrom ggpubr ggarrange
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off tiff
#' @importFrom foreach foreach %dopar% %do%
MK_ProtConnMult <- function(nodes, regions, thintersect = NULL,
                            attribute = "Intersected area", area_unit = "ha",
                            res_attribute = 150,
                            distance = list(type= "centroid", resistance = NULL),
                            distance_thresholds, probability,
                            transboundary = NULL,  keep = 1,
                            CI = "all",
                            plot = FALSE,
                            write = NULL, intern = TRUE,
                            parallel = FALSE){
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (missing(regions)) {
    stop("error missing file of region")
  } else {
    if (is.numeric(regions) | is.character(regions)) {
      stop("error missing file of regions")
    }
  }

  if (is.null(probability) | !is.numeric(probability)) {
    stop("error missing probability")
  }

  if (is.null(distance_thresholds)) {
    stop("error missing numeric distance threshold(s)")
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(nrow(regions)<=1){
    stop("one region, for better results please use MK_ProtConn()")
  }

  #Remove error zm
  if(class(regions)[1] == "SpatialPolygonsDataFrame") {
    regions <- st_as_sf(regions)
  }
  regions <- st_zm(regions)
  regions$ID_Temp <- 1:nrow(regions)

  if(class(nodes)[1] == "SpatialPolygonsDataFrame") {
    nodes <- st_as_sf(nodes)
  }
  nodes <- st_zm(nodes)

  if(isFALSE(parallel)){
    pb <- progress_estimated(length(regions$ID_Temp), 0)
    x=NULL
    protconn_result <- map(1:length(regions$ID_Temp), function(x){
      Ecoreg_sel <- regions[regions$ID_Temp == unique(regions$ID_Temp)[x],]
      if (isTRUE(intern)){
        pb$tick()$print()
      }
      protconn <- tryCatch(MK_ProtConn(nodes = nodes,
                                       region = Ecoreg_sel,
                                       thintersect = thintersect,
                                       attribute = attribute,
                                       area_unit = area_unit,
                                       distance = distance,
                                       res_attribute = res_attribute,
                                       transboundary = transboundary,
                                       keep = keep,
                                       distance_thresholds = distance_thresholds,
                                       probability = probability,
                                       plot = FALSE, intern = FALSE),  error = function(err)err)


      if (inherits(protconn, "error")){
        stop(protconn)
      }

      if(length(distance_thresholds) > 1){
        Ecoreglist <- tryCatch(map(protconn, function(x){
          n <- as.vector(x[[1]][[3]])
          x2 <- t(x[[1]][[4]]) %>% as.data.frame()
          colnames(x2) <- n
          return(x2)}), error = function(err)err)
      } else {
        n <- as.vector(protconn[[3]])
        Ecoreglist <- t(protconn[[4]]) %>% as.data.frame()
        colnames(Ecoreglist) <- n
      }


      if (inherits(Ecoreglist, "error")){
        stop(Ecoreglist)
      }

      if(length(distance_thresholds) > 1){
      ECAlist <- tryCatch(map(protconn, function(x){
        n <- as.vector(x[[1]][[1]])
        x2 <- t(x[[1]][[2]]) %>% as.data.frame()
        colnames(x2) <- n
        return(x2)}), error = function(err)err)
      } else {
        n <- as.vector(protconn[[1]])
        ECAlist <- t(protconn[[2]]) %>% as.data.frame()
        colnames(ECAlist) <- n
      }

      if (inherits(ECAlist, "error")){
        stop(ECAlist)
      }

      result <- list(Ecoreglist, ECAlist)
      return(result)})
  } else {
    plan(strategy = multiprocess)
    protconn_result <- tryCatch(future_map(as.list(1:length(regions$ID_Temp)), function(x){
      Ecoreg_sel <- regions[regions$ID_Temp == unique(regions$ID_Temp)[x],]
      protconn <- tryCatch(MK_ProtConn(nodes = nodes,
                                       region = Ecoreg_sel,
                                       thintersect = thintersect,
                                       attribute = attribute,
                                       area_unit = area_unit,
                                       res_attribute = res_attribute,
                                       distance = distance,
                                       transboundary = transboundary,
                                       keep = keep,
                                       distance_thresholds = distance_thresholds,
                                       probability = probability,
                                       plot = FALSE, intern = FALSE),  error = function(err)err)

      if (inherits(protconn, "error")){
        stop(protconn)
      }

      if(length(distance_thresholds) > 1){
        Ecoreglist <- tryCatch(map(protconn, function(x){
          n <- as.vector(x[[1]][[3]])
          x2 <- t(x[[1]][[4]]) %>% as.data.frame()
          colnames(x2) <- n
          return(x2)}), error = function(err)err)
      } else {
        n <- as.vector(protconn[[3]])
        Ecoreglist <- t(protconn[[4]]) %>% as.data.frame()
        colnames(Ecoreglist) <- n
      }


      if (inherits(Ecoreglist, "error")){
        stop(Ecoreglist)
      }

      if(length(distance_thresholds) > 1){
        ECAlist <- tryCatch(map(protconn, function(x){
          n <- as.vector(x[[1]][[1]])
          x2 <- t(x[[1]][[2]]) %>% as.data.frame()
          colnames(x2) <- n
          return(x2)}), error = function(err)err)
      } else {
        n <- as.vector(protconn[[1]])
        ECAlist <- t(protconn[[2]]) %>% as.data.frame()
        colnames(ECAlist) <- n
      }


      if (inherits(ECAlist, "error")){
        stop(ECAlist)
      }

      result <- list(Ecoreglist, ECAlist)
      return(result)}, .progress = intern), error = function(err)err)
  }

  Ecoreglist <- lapply(protconn_result, function(x){x[[1]]})

  ECAlist <- lapply(protconn_result, function(x){x[[2]]})

  i=NULL
  j=NULL
  data_1 <- tryCatch(map(1:length(distance_thresholds), function(i){

    if(length(distance_thresholds)>1){
      data_2 <- map(1:length(regions$ID_Temp), function(j){
        Ecoreglist[[j]][[i]]
      })
    } else {
      data_2 <- Ecoreglist
    }

    ProtConnEcor <- do.call(rbind, data_2)
    rownames(ProtConnEcor) <- NULL
    ProtConnEcor$ID_Temp <- regions$ID_Temp

    if(length(distance_thresholds)>1){
      ECA_2 <- map(1:length(regions$ID_Temp), function(j) {
        ECAlist[[j]][[i]]
      })
      } else {
        ECA_2 <- ECAlist
      }

    ECA_3 <- do.call(rbind, ECA_2)
    rownames(ECA_3) <- NULL
    ECA_3 <- ECA_3[,c(1:2)]
    ECA_3$ID_Temp <- regions$ID_Temp

    ProtConnEcor <- base::merge(ECA_3, ProtConnEcor, by = "ID_Temp")
    regions.2 <- base::merge(regions, ProtConnEcor, by = "ID_Temp")
    regions.2$ID_Temp <- NULL
    results <- list()
    results[[2]] <- regions.2
    DataProtconn <- regions.2
    st_geometry(DataProtconn) <- NULL

    #
    DataProtconn_2 <-  tryCatch(ProtConnStat(DataProtconn, ci = CI, nr = 1000), error = function(err)err)
    if(inherits(DataProtconn_2, "error")){
      DataProtconn_2 <-  tryCatch(ProtConnStat(DataProtconn, ci = NULL), error = function(err)err)
      if(inherits(DataProtconn_2, "error")){
        DataProtconn_2 <- DataProtconn
      }
    }
    DataProtconn_2$Indicator <- rownames(DataProtconn_2)
    DataProtconn_2 <- DataProtconn_2[moveme(names(DataProtconn_2), "Indicator first")]
    rownames(DataProtconn_2) <- NULL
    names(DataProtconn_2)[1:2] <- c("ProtConn indicator", "Values(%)")
    DataProtconn_2[c(2:ncol(DataProtconn_2))] <- round(DataProtconn_2[c(2:ncol(DataProtconn_2))], 3)

    if(!is.null(write)){
      write.csv(DataProtconn, paste0(write, "SummaryStats_", distance_thresholds[i], ".csv"))
      }

    if(ncol(DataProtconn_2) > 2){
      DataProtconn <- formattable(DataProtconn_2[3:nrow(DataProtconn_2),], align = c("l", rep("r", NCOL(DataProtconn_2) - 1)),
                                  list(`ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                       `Values(%)` = color_tile("white", "#F88B13"),
                                       area(col = 3:4) ~ color_tile("white", "#CE5D9B")))
    } else {
      DataProtconn <- formattable(DataProtconn_2[3:nrow(DataProtconn_2),], align = c("l", rep("r", NCOL(DataProtconn_2) - 1)),
                                  list(`ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                       `Values(%)` = color_tile("white", "#F88B13")))
    }

    results[[1]] <- DataProtconn
    names(results) <- c(paste0("ProtConn_", distance_thresholds[i]),
                        paste0("ProtConn_overall_", distance_thresholds[i]))

    if(isTRUE(plot)){
        dacc <- as.data.frame(DataProtconn[c(1:3), c(2:3)])
        dacc$name <- c("Unprotected", "Protected", "Protected connected")
        dacc$name <- factor(dacc$name, levels = c("Unprotected", "Protected", "Protected connected"))
        dacc$col <- c("#C34D51", "#53A768", "#4C72AF")
        names(dacc)[1] <- "Values"
        dacc <- dacc[which(dacc$Values > 0), ]
        if(nrow(dacc) > 1){
          dacc$max <- dacc$Values+dacc$SD
          dacc$min <- dacc$Values-dacc$SD
          dacc[5:6] <- apply(dacc[5:6], 2, function(x){
            x[which(x >100)]<- 100
            x[which(x <0)]<- 0
            return(x)})

          plot_protconn1 <- ggplot(dacc, aes(x = dacc$name, y = dacc$Values, fill = dacc$name)) +
            geom_bar(position = position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
            geom_errorbar(position = position_dodge(), width = 0.2, aes(ymin = min, ymax = max)) +
            labs(title = paste0("ProtConn Indicators: ", distance_thresholds[i]), x = "", y = "Percentage (%)", size = rel(1.2)) +
            theme_bw()  +
            theme(plot.title = element_text(color = "#252525", size = rel(1.4), hjust = 0.5, face = "bold"),
                  axis.title= element_text(color = "#252525", size = rel(1.2)),
                  legend.title= element_blank(),
                  legend.text = element_text(colour = "#252525", size = rel(1.2)),
                  axis.text= element_text(colour = "#525252", size = rel(1)))+
            scale_fill_manual(values = dacc$col) +
            geom_hline(aes(yintercept = 17, linetype = "Aichi Target (17%)"),
                       colour = 'black', size = 1.2) +
            scale_linetype_manual(name = " Aichi Target", values = c(2, 2),
                                  guide = guide_legend(override.aes = list(color = c("black"), size = 0.8)))
          plots <- list(plot_protconn1)
        }

        dacc2 <- as.data.frame(DataProtconn[c(7:10), c(1:3)])
        dacc2$name <- c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]")
        dacc2$Indicator <- NULL
        dacc2$name <- factor(dacc2$name, levels = c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]"))
        dacc2$col <- c("#253494", "#2c7fb8", "#41b6c4", "#7fcdbb")
        names(dacc2)[2] <- "Values"
        dacc2 <- dacc2[which(dacc2$Values > 0), ]
        if(nrow(dacc2) > 1){
          dacc2$max <- dacc2$Values+dacc2$SD
          dacc2$min <- dacc2$Values-dacc2$SD
          dacc2[6:7] <- apply(dacc2[6:7], 2, function(x){
            x[which(x >100)]<- 100
            x[which(x <0)]<- 0
            return(x)})
          plot_protconn2 <- ggplot(dacc2, aes(x = dacc2$name, y = dacc2$Values, fill = dacc2$name)) +
            geom_bar(position = position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
            geom_errorbar(position = position_dodge(), width = 0.2, aes(ymin = min, ymax = max)) +
            labs(title = "Protected connected fraction", x = "", y = "Percentage (%)", size = rel(1.2)) +
            theme_bw()  +
            theme(plot.title = element_text(color = "#252525", size = rel(1.4), hjust = 0.45, face = "bold"),
                  axis.title= element_text(color = "#252525", size = rel(1.2)),
                  plot.margin = margin(0, 5.3, 0, 0.3, "cm"),
                  legend.title= element_blank(),
                  legend.text = element_text(colour = "#252525", size = rel(1.2)),
                  axis.text= element_text(colour = "#525252", size = rel(1)))+
            scale_fill_manual(values = dacc2$col)
          plots <- list(plot_protconn1, plot_protconn2)
        }
        plots <- compact(plots)

        if(length(plots) == 2){
          figure <- ggarrange(plots[[1]], plots[[2]],
                              ncol = 1, nrow = 2)
        } else if (length(plots) == 1) {
          figure <- plots[[1]]
        } else {
          figure <- "There are insufficient data to plot"
        }

        results[[3]] <- figure
        names(results)[3] <- "ProtConn Plot"
      }

    if(!is.null(write)){
        regions.2[is.na(regions.2)] <- 0
        nn <- c("EC", "PC", "Unprot", "Prot", "ProtConn", "ProtUnconn", "RelConn", "Design", "Bound", "P_Prot",
                "P_Trans", "P_Unprot", "P_Within", "P_Contig", "P_WithinL", "P_ContigL", "P_UnprotL", "P_TransL")
        names(regions.2)[which(names(regions.2)=="EC(PC)"):(ncol(regions.2)-1)] <- nn

        write_sf(regions.2, paste0(write, "ProtConn_", distance_thresholds[i], ".shp"), delete_layer = T)
        if(isTRUE(plot)){
          if(!is.character(figure)){
            tiff(paste0(write, "ProtConn_plot_d", distance_thresholds[i], ".tif"), width = 806, height = 641)
            print(figure)
            dev.off()
          }
        }
      }

    return(results)
    }), error = function(err)err)

  if (inherits(data_1, "error")){
    stop(data_1)
  }
  names(data_1) <- paste0("ProtConn_", distance_thresholds)
  return(data_1)
  }
