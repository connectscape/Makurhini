#' Multiple Protected Connected (ProtConn)
#'
#' Estimate Protected Connected (ProtConn) indicator and fractions for multiple regions.
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param regions object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param area_unit character. Attribute area units. You can set an area unit, "Makurhini::unit_covert()" compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to hectares "m2".
#' @param distance list. See \link[Makurhini]{distancefile}. Example, list(type= "centroid", resistance = NULL).
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_threshold = c(30000, 50000); sequence distances: distance_threshold = seq(10000,100000, 10000).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. numeric.
#' If probability = NULL, then it will be the inverse of the mean dispersal distance for the species (1/α; Hanski and Ovaskainen 2000).
#' @param transboundary numeric. Buffer to select polygons in a second round, their attribute value = 0, see Saura et al. 2017.
#' @param transboundary_type character. Two options: "nodes" or "region". If it is "nodes" the transboundary is built from the limits of the nodes present in the region (default),
#' if "region" is selected the transboundary is built from the limits of the region.
#' @param protconn_bound logical. If TRUE then the fractions ProtUnConn[design] and ProtConn[bound] will be estimated.
#' @param geom_simplify logical. Slightly simplify the region and nodes geometries.
#' @param delta logical. Estimate the contribution of each node to the ProtConn value in the regions.
#' @param CI character. A character vector representing the type of confidence intervals that will be estimated. The value should be any subset of the values c("norm","basic", "stud", "perc", "bca") or "all" which will compute all five types of intervals (see, \link[boot]{boot.ci})
#' @param plot logical. Plot the main ProtConn indicators and fractions with their standard deviation, default = FALSE.
#' @param write character. Output folder including the output file name without extension, e.g., "C:/ProtConn/Protfiles".
#' @param intern logical. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @param parallel numeric. Specify the number of cores to use for parallel processing, default = NULL.
#' Parallelize the function using furrr package and multiprocess plan.
#' @return
#' Table with the following ProtConn values: ECA, Prot, ProtConn, ProtUnconn, RelConn, ProtUnConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within],
#'  ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot], ProtConn[Within][land], ProtConn[Contig][land],
#'  ProtConn[Unprot][land], ProtConn[Trans][land] \cr
#' \cr
#' *Each indicator is accompanied by six dispersion statistics: standard deviation, standard error and four confidence intervals obtained with a bootstrap-type resampling (Carpenter & Bithell 2000; see, \strong{ProtConnStat()})\cr
#' \cr
#' *If plot is not TRUE a list is returned with the ProtConn table and a plots.
#' @references Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S., Bertzky, B., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2018). Protected area connectivity: Shortfalls in global targets and country-level priorities. Biological Conservation, 219(October 2017), 53–67.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(raster)
#' data("Protected_areas", package = "Makurhini")
#' plot(Protected_areas, col="green")
#'
#' data("regions", package = "Makurhini")
#' plot(regions, col=c("blue", "red", "green"))
#'
#' test <- MK_ProtConnMult(nodes = Protected_areas,
#'                         regions = regions,
#'                         area_unit = "ha",
#'                         distance = list(type= "centroid"),
#'                         distance_thresholds = c(10000, 50000),
#'                         probability = 0.5, transboundary = 50000,
#'                         plot = TRUE, write = NULL,
#'                         parallel = NULL, intern = TRUE)
#' test
#' }
#' @importFrom magrittr %>%
#' @importFrom sf st_as_sf st_zm st_geometry st_geometry<- write_sf
#' @importFrom progressr handlers handler_pbcol progressor
#' @importFrom crayon bgWhite white bgCyan
#' @importFrom future multiprocess plan availableCores
#' @importFrom furrr future_map future_map_dfr
#' @importFrom formattable color_tile area style formattable formatter
#' @importFrom utils write.csv
#' @importFrom ggplot2 ggplot geom_bar aes position_dodge labs rel theme_bw theme element_blank element_text scale_fill_manual geom_hline scale_linetype_manual guide_legend margin geom_errorbar
#' @importFrom purrr compact map_df
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices dev.off tiff
MK_ProtConnMult <- function(nodes, regions,
                            area_unit = "m2",
                            distance = list(type= "centroid", resistance = NULL),
                            distance_thresholds, probability,
                            transboundary = NULL,
                            transboundary_type = "nodes",
                            protconn_bound = FALSE,
                            geom_simplify = FALSE,
                            delta = FALSE,
                            CI = "all",
                            plot = FALSE,
                            write = NULL, intern = TRUE,
                            parallel = NULL){
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

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(nrow(regions)<=1){
    stop("one region, for better results please use MK_ProtConn()")
  } else {
    regions$ID_Temp <- 1:nrow(regions)
  }

  regions <- TopoClean(regions, xsimplify = geom_simplify)
  loop <- 1:nrow(regions)

  if (isTRUE(intern)){
    handlers(global = T)
    handlers(handler_pbcol(complete = function(s) crayon::bgYellow(crayon::white(s)),
                           incomplete = function(s) crayon::bgWhite(crayon::black(s)),
                           intrusiveness = 2))
    pb <- progressor(along = loop)
    message("Step 1. ProtConn estimation")
  }


  if(is.null(parallel)){
    protconn_result <- tryCatch(lapply(loop, function(x){
      Ecoreg_sel <- regions[x,]

      if (isTRUE(intern)){
        pb()
      }

      protconn <- tryCatch(MK_ProtConn(nodes = nodes,
                                       region = Ecoreg_sel,
                                       area_unit = area_unit,
                                       distance = distance,
                                       transboundary = transboundary,
                                       transboundary_type = transboundary_type,
                                       protconn_bound = protconn_bound,
                                       distance_thresholds = distance_thresholds,
                                       probability = probability,
                                       geom_simplify = geom_simplify,
                                       delta = delta,
                                       plot = FALSE,
                                       intern = FALSE),  error = function(err)err)


      if (inherits(protconn, "error")){
        stop(protconn)
      }
      return(protconn)
    }), error = function(err)err)

    if (inherits(protconn_result, "error")){
      stop("Error, please review your input shapefiles")
    }

    if(length(distance_thresholds)>1){
      Ecoreglist <- tryCatch(lapply(1:length(distance_thresholds), function(x){
        x.1 <- map_df(protconn_result, function(y){
          y.1 <- y[[x]]
          y.2 <- y.1[[3]]
          y.3 <- t(y.1[[4]]) %>% as.data.frame()
          colnames(y.3) <- y.2
          return(y.3)
        })
        return(x.1)
      }), error = function(err)err)

      ECAlist <- tryCatch(lapply(1:length(distance_thresholds), function(x){
        x.1 <- map_df(protconn_result, function(y){
          y.1 <- y[[x]]
          y.2 <- as.vector(y.1[[1]][1:2])
          y.3 <- t(as.numeric(y.1[[2]][1:2])) %>% as.data.frame()
          colnames(y.3) <- y.2
          return(y.3)
        })
        return(x.1)}), error = function(err)err)

      if (inherits(Ecoreglist, "error")|inherits(ECAlist, "error")){
        stop("Errorc in table construction")
      }
    } else {
      Ecoreglist <- tryCatch(map_df(protconn_result, function(x){
        x.1 <- x[[3]]
        x.2 <- t(x[[4]]) %>% as.data.frame()
        colnames(x.2) <- x.1
        return(x.2)
      }), error = function(err)err)

      ECAlist <- tryCatch(map_df(protconn_result, function(x){
        x.1 <- as.vector(x[[1]][1:2])
        x.2 <- t(x[[2]][1:2]) %>% as.data.frame()
        colnames(x.2) <- x.1
        return(x.2)
      }), error = function(err)err)

      if (inherits(Ecoreglist, "error")|inherits(ECAlist, "error")){
        stop("Error in table construction")
      } else {
        Ecoreglist <- list(Ecoreglist)
        ECAlist <- list(ECAlist)
      }
    }

  } else {
    works <- as.numeric(availableCores())-1
    works <- if(parallel > works){works}else{parallel}
    plan(strategy = multiprocess, gc = TRUE, workers = works)
    protconn_result <- tryCatch(future_map(loop, function(x){
      if (isTRUE(intern)){
        pb()
      }
      Ecoreg_sel <- regions[regions$ID_Temp == unique(regions$ID_Temp)[x],]
      protconn <- tryCatch(MK_ProtConn(nodes = nodes,
                                       region = Ecoreg_sel,
                                       area_unit = area_unit,
                                       distance = distance,
                                       transboundary = transboundary,
                                       transboundary_type = transboundary_type,
                                       distance_thresholds = distance_thresholds,
                                       probability = probability,
                                       geom_simplify = geom_simplify,
                                       delta = delta,
                                       plot = FALSE,
                                       intern = FALSE),  error = function(err)err)


      if (inherits(protconn, "error")){
        stop(protconn)
      }

      return(protconn)}), error = function(err)err)

    if(inherits(protconn_result, "error")){
      close_multiprocess(works)
      stop("Error, please review your input shapefiles")
    }

    if(length(distance_thresholds)>1){
      Ecoreglist <- tryCatch(future_map(1:length(distance_thresholds), function(x){
        x.1 <- future_map_dfr(protconn_result, function(y){
          y.1 <- y[[x]]
          y.2 <- y.1[[3]]
          y.3 <- t(y.1[[4]]) %>% as.data.frame()
          colnames(y.3) <- y.2
          return(y.3)
        })
        return(x.1)
      }), error = function(err)err)

      ECAlist <- tryCatch(future_map(1:length(distance_thresholds), function(x){
        x.1 <- future_map_dfr(protconn_result, function(y){
          y.1 <- y[[x]]
          y.2 <- as.vector(y.1[[1]][1:2])
          y.3 <- t(as.numeric(y.1[[2]][1:2])) %>% as.data.frame()
          colnames(y.3) <- y.2
          return(y.3)
        })
        return(x.1)}), error = function(err)err)

      if (inherits(Ecoreglist, "error")|inherits(ECAlist, "error")){
        close_multiprocess(works)
        stop("Error in table construction")
      }
    } else {
      Ecoreglist <- tryCatch(future_map_dfr(protconn_result, function(x){
        x.1 <- x[[3]]
        x.2 <- t(x[[4]]) %>% as.data.frame()
        colnames(x.2) <- x.1
        return(x.2)
      }), error = function(err)err)

      ECAlist <- tryCatch(future_map_dfr(protconn_result, function(x){
        x.1 <- as.vector(x[[1]][1:2])
        x.2 <- t(x[[2]][1:2]) %>% as.data.frame()
        colnames(x.2) <- x.1
        return(x.2)
      }), error = function(err)err)

      if (inherits(Ecoreglist, "error")|inherits(ECAlist, "error")){
        close_multiprocess(works)
        stop("Error in table construction")
      } else {
        Ecoreglist <- list(Ecoreglist)
        ECAlist <- list(ECAlist)
      }
    }

    close_multiprocess(works)
  }

  if (isTRUE(intern)){
    message("Step 2. Estimating statistics")
  }

  protconn_result2 <- tryCatch(lapply(1:length(distance_thresholds), function(x){
    x.1 <- Ecoreglist[[x]]
    rownames(x.1) <- NULL
    x.1$ID_Temp <- regions$ID_Temp

    x.2 <- ECAlist[[x]]
    rownames(x.2) <- NULL
    x.2$ID_Temp <- regions$ID_Temp

    x.3 <- base::merge(x.2, x.1, by = "ID_Temp")
    x.3 <- base::merge(regions, x.3, by = "ID_Temp")
    x.3$ID_Temp <- NULL

    results <- list()
    results[[2]] <- x.3

    if(!is.null(write)){
      x.0 <- x.3
      x.0[is.na(x.0)] <- 0

      if(isTRUE(protconn_bound)){
        nn <- c("EC", "PC", "Prot", "Unprot", "ProtConn", "ProtUnconn",
                "Design", "Bound", "RelConn", "P_Prot",
                "P_Trans", "P_Unprot", "P_Within", "P_Contig", "P_WithinL", "P_ContigL", "P_UnprotL", "P_TransL")
      } else {
        nn <- c("EC", "PC", "Prot", "Unprot", "ProtConn", "ProtUnconn",
                "RelConn", "P_Prot", "P_Trans", "P_Unprot", "P_Within", "P_Contig", "P_WithinL", "P_ContigL",
                "P_UnprotL", "P_TransL")
      }

      names(x.0)[which(names(x.0)=="EC(PC)"):(ncol(x.0)-1)] <- nn

      write_sf(x.0, paste0(write, "ProtConn_", distance_thresholds[x], ".shp"), delete_layer = T)
    }


    st_geometry(x.3) <- NULL

    if(is.null(CI)){
      x.4 <-  tryCatch(ProtConnStat(x.3, ci = NULL), error = function(err)err)
      if(inherits(x.4, "error")){
        x.4 <- x.3
      }
    } else {
      x.4 <-  tryCatch(ProtConnStat(x.3, ci = CI, nr = 1000), error = function(err)err)

      if(inherits(x.4, "error")){
        x.4 <-  tryCatch(ProtConnStat(x.3, ci = NULL), error = function(err)err)

        if(inherits(x.4, "error")){
          x.4 <- x.3
        }
      }
    }

    x.4$Indicator <- rownames(x.4)
    x.4 <- x.4[moveme(names(x.4), "Indicator first")]
    rownames(x.4) <- NULL
    names(x.4)[1:2] <- c("ProtConn indicator", "Values (%)")
    x.4[c(2:ncol(x.4))] <- round(x.4[c(2:ncol(x.4))], 3)

    if(!is.null(write)){
      write.csv(x.4, paste0(write, "SummaryStats_", distance_thresholds[x], ".csv"), row.names = F)
    }

    if(ncol(x.4) > 2){
      row.names(x.4) <- NULL
      DataProtconn <- formattable(x.4[3:nrow(x.4),], align = c("l", rep("r", NCOL(x.4) - 1)),
                                  list(`ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                       `Values (%)` = color_tile("white", "#F88B13"),
                                       area(col = 3:4) ~ color_tile("white", "#CE5D9B")))
    } else {
      row.names(x.4) <- NULL
      DataProtconn <- formattable(x.4[3:nrow(x.4),], align = c("l", rep("r", NCOL(x.4) - 1)),
                                  list(`ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                       `Values (%)` = color_tile("white", "#F88B13")))
    }

    results[[1]] <- DataProtconn
    names(results) <- c(paste0("ProtConn_overall", distance_thresholds[x]),
                        paste0("ProtConn_", distance_thresholds[x]))
    return(results)}), error = function(err)err)

  if (inherits(protconn_result2, "error")){
    stop("Error in CI estimation")
  }


  if(isTRUE(plot)){
    if(isTRUE("ggplot2" %in% rownames(installed.packages())) &
       isTRUE("ggpubr" %in% rownames(installed.packages()))){
      protconn_result2 <- tryCatch(map(1:length(distance_thresholds), function(x){
        DataProtconn <- protconn_result2[[x]][[1]] %>% as.data.frame()
        dacc <- DataProtconn[c(2,1,3), 2:3]
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
            ggplot2::geom_errorbar(position = position_dodge(), width = 0.2, aes(ymin = min, ymax = max)) +
            labs(title = paste0("ProtConn Indicators: ", distance_thresholds[x]), x = "", y = "Percentage (%)", size = rel(1.2)) +
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

        dacc2 <- DataProtconn[c(7:10), c(1:3)]
        dacc2$name <- c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]")
        dacc2$Indicator <- NULL
        dacc2$name <- factor(dacc2$name, levels = c("ProtConn[Within]", "ProtConn[Contig]", "ProtConn[Unprot]", "ProtConn[Trans]"))
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
            ggplot2::geom_errorbar(position = position_dodge(), width = 0.2, aes(ymin = min, ymax = max)) +
            labs(title = "Protected connected fractions", x = "", y = "Percentage (%)", size = rel(1.2)) +
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
          figure <- ggarrange(plots[[1]], plots[[2]], ncol = 1, nrow = 2)
        } else if (length(plots) == 1) {
          figure <- plots[[1]]
        } else {
          figure <- "There are insufficient data to plot"
        }

        if(!is.null(write)){
          if(!is.character(figure)){
            tiff(paste0(write, "ProtConn_plot_d", distance_thresholds[x], ".tif"), width = 806, height = 641)
            print(figure)
            dev.off()
          }
        }

        protconn_result2[[x]][[3]] <- figure
        names(protconn_result2[[x]])[3] <- "ProtConn Plot"
        return(protconn_result2[[x]])
      }), error = function(err)err)
      if (inherits(protconn_result2, "error")){
        stop("Error plotting ProtConn")
      }
    } else {
      message("To make the plots you need to install the packages ggplot2 and ggpubr")
    }
  }
  names(protconn_result2) <- paste0("ProtConn_", distance_thresholds)
  return(protconn_result2)
}

