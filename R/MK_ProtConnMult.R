#' Multiple Protected Connected (ProtConn)
#'
#' Estimate Protected Connected (ProtConn) indicator and fractions for multiple regions.
#' @param nodes object of class \code{sf, sfc, sfg, spatialPolygonsDataFrame}. Spatial data of vector type that normally contains the spatial limits of protected areas. It must be in a projected coordinate system.
#' @param regions object of class \code{sf, sfc, sfg, spatialPolygonsDataFrame}. Polygon delimiting the \bold{regions} or \bold{study areas}. It must be in a projected coordinate system.
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
#' @param delta \code{logical}. Estimate the contribution of each node to the ProtConn value in each region.
#' @param CI \code{character}. A character vector representing the type of confidence intervals that will be estimated. The value should be any subset of the values \code{c("norm","basic", "stud", "perc", "bca")} or \code{"all"} which will compute all five types of intervals (see, \link[boot]{boot.ci})
#' @param plot \code{logical}. Plot the main ProtConn indicators and fractions with their standard deviation, default = \code{FALSE}.
#' @param parallel \code{numeric}. Specify the number of cores to use for parallel processing, default = NULL. Parallelize the function using furrr package and multiprocess plan.
#' @param write \code{character}. Output folder including the output file name without extension, e.g., \code{"C:/ProtConn/Protfiles"}.
#' @param intern \code{logical}. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @return
#' For each region:\cr
#' - Table with the following ProtConn values: \code{ECA, Prot, ProtConn, ProtUnconn, RelConn, ProtUnConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within],
#'  ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot], ProtConn[Within][land], ProtConn[Contig][land],
#'  ProtConn[Unprot][land], ProtConn[Trans][land]} \cr
#' \cr
#' - If plot \bold{is not NULL} a \code{list} is returned with the ProtConn table and a plots.
#' - If \bold{delta} is \code{TRUE} then it returns an sf class object with the importance value (contribution to ProtConn) for each node in the region.
#' @references Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S., Bertzky, B., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2018). Protected area connectivity: Shortfalls in global targets and country-level priorities. Biological Conservation, 219(October 2017), 53–67.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(sf)
#' data("Ecoregions", package = "Makurhini")#'
#' #For this example, we select the first three columns and the first 10 of the ecoregions
#' Ecoregions <- Ecoregions[1:15,1:3]
#' plot(st_geometry(Ecoregions), col = "#7E6A9F")
#' load(system.file("extdata", "Protected_areas.rda",
#'                 package = "Makurhini", mustWork = TRUE))
#' #plot(st_geometry(Protected_areas), col="green", add = TRUE) #It may take time to plot all PAs.
#'
#' test <- MK_ProtConnMult(nodes = Protected_areas,
#'                         regions = Ecoregions,
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
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map future_map_dfr
#' @importFrom formattable color_tile area style formattable formatter
#' @importFrom utils write.csv
#' @importFrom ggplot2 ggplot geom_bar aes position_dodge labs rel theme_bw theme element_blank element_text scale_fill_manual geom_hline scale_linetype_manual guide_legend margin geom_errorbar
#' @importFrom purrr compact map_df
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices dev.off tiff
#' @importFrom utils installed.packages
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
                            write = NULL,
                            parallel = NULL,
                            intern = TRUE){
  if(isTRUE(intern)){
    message("Step 1. Reviewing parameters")
  }
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

  if(isFALSE(parallel)){
    parallel <- NULL
  }

  if(isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
  }

  if(nrow(regions)<=1){
    stop("one region, for better results please use MK_ProtConn()")
  } else {
    regions$ID_Temp <- 1:nrow(regions)
  }

  regions <- TopoClean(regions, xsimplify = geom_simplify)
  if (isTRUE(intern)){
    message("Step 2. Processing ProtConn metric. Progress estimated:")
    pb <- txtProgressBar(0, nrow(regions), style = 3)
  }


  if(is.null(parallel)){
    protconn_result <- tryCatch(lapply(1:nrow(regions), function(x){
      Ecoreg_sel <- regions[x,]

      if (isTRUE(intern)){
        setTxtProgressBar(pb, x)
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
    works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}
    if(.Platform$OS.type == "unix") {
      strat <- future::multicore
    } else {
      strat <- future::multisession
    }
    plan(strategy = strat, gc = TRUE, workers = works)
    protconn_result <- tryCatch(future_map(1:nrow(regions), function(x){
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

      return(protconn)}, .progress = intern), error = function(err)err)

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
                                       formattable::area(col = 3:4) ~ color_tile("white", "#CE5D9B")))
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
        dacc$col <- c("#fc8d62", "#66c2a5", "#8da0cb")
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
            geom_hline(aes(yintercept = 17, linetype = "Aichi Target (17%)"), colour = 'black', size = 1.2) +
            geom_hline(aes(yintercept = 30, linetype = "Kunming-Montreal (30%)"), colour = 'red', size = 1.2)+
            scale_linetype_manual(name = " Aichi Target", values = c(2, 2),
                                  guide = guide_legend(override.aes = list(color = c("black", 'red'), size = 0.8)))
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

