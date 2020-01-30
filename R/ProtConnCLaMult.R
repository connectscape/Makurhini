#' Multiple Protected Connected (ProtConn)
#'
#' Use the CONEFOR command line to estimate Protected Connected (ProtConn) indexes for multiple ecoregions
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. Protected areas (PAs) shapefile.
#' @param regions object of class sf, sfc, sfg or SpatialPolygons. Ecoregions shapefile.
#' @param thintersect numeric.Threshold of intersection in percentage allowed to select or not a target geometry. Default = 90, if intersection >=90 percentage, the geometry will be selected.
#' @param attribute character. Select the nodes attribute: "Area" = Complete Protected areas; "Intersected area" = Intersected Protected areas; or another specific column name with the nodes attribute, ideally this attribute mus be an area-weighted index, otherwise the interpretation of the protconn index may change.
#' @param res_attribute numeric. If the attribute is no equal to "Area" or "Intersected area" then nodes will be converted to raster to extract values in one  process step, you can set the raster resolution, default = 150.
#' @param distance list.See distancefile(). E.g.: list(type = "centroid", resistance = NULL, geometry_out = NULL).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' @param distance_threshold numeric. Distance or distances thresholds to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_threshold = c(30000, 50000); sequence distances: distance_threshold = seq(10000,100000, 10000).
#' @param transboundary numeric. Buffer to select polygones in a second round, their attribute value = 0, see Saura et al. 2017.
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
ProtConnCLaMult <- function(nodes, regions, thintersect = NULL,
                            attribute = "Intersected area", res_attribute = 150,
                            distance = list(type= "centroid", resistance = NULL, tolerance = NULL, geometry_out = NULL, ...),
                            distance_thresholds, probability,
                            transboundary = NULL,
                            CI = "all",
                            plot = FALSE,
                            write = NULL, intern = TRUE,
                            parallel = FALSE){
  #Remove error zm
  regions <- sf::st_as_sf(regions) %>% sf::st_zm(.)
  #Nuevo ID a regiones
  regions$ID_Temp <- 1:nrow(regions)
  nodes <- sf::st_zm(nodes)
  ttt.2 <- getwd()

  #
  if(isFALSE(parallel)){
    pb <- dplyr::progress_estimated(length(regions$ID_Temp), 0)

    protconn_result <- tryCatch(purrr::map(as.list(1:length(regions$ID_Temp)), function(x){
      Ecoreg_sel <- regions[regions$ID_Temp == unique(regions$ID_Temp)[x],]
      if (isTRUE(intern)){
        pb$tick()$print()
        }
    protconn <- tryCatch(ProtConnCLa(nodes = nodes,
                                     region = Ecoreg_sel,
                                     thintersect = thintersect,
                                     attribute = attribute,
                                     distance = distance,
                                     transboundary = transboundary,
                                     distance_thresholds = distance_thresholds,
                                     probability = probability,
                                     plot = FALSE, intern = FALSE),  error = function(err)err)

    if (inherits(protconn, "error")){
      stop(protconn)
    }

    Ecoreglist <- tryCatch(purrr::map(protconn, function(x){
      n <- as.vector(x[[1]][[3]])
      x2 <- t(x[[1]][[4]]) %>% as.data.frame()
      colnames(x2) <- n
      return(x2)}), error = function(err)err)

    if (inherits( Ecoreglist, "error")){
      Ecoreglist <- purrr::map(protconn, function(x){
        n <- as.vector(x[[1]][[3]])
        x2 <- t(x[[1]][[4]]) %>% as.data.frame()
        colnames(x2) <- n
        return(x2)})
    }

    ECAlist <- tryCatch(purrr::map(protconn, function(x){
      n <- as.vector(x[[1]][[1]])
      x2 <- t(x[[1]][[2]]) %>% as.data.frame()
      colnames(x2) <- n
      return(x2)}), error = function(err)err)

    if (inherits(ECAlist, "error")){
      ECAlist <- purrr::map(protconn, function(x){
        n <- as.vector(x[[1]][[1]])
        x2 <- t(x[[1]][[2]]) %>% as.data.frame()
        colnames(x2) <- n
        return(x2)})
    }
    result <- list(Ecoreglist, ECAlist)
    return(result)}), error = function(err)err)
  } else {
    future::plan(strategy = future::multiprocess)
    protconn_result <- tryCatch(furrr::future_map(as.list(1:length(regions$ID_Temp)), function(x){
      Ecoreg_sel <- regions[regions$ID_Temp == unique(regions$ID_Temp)[x],]
      protconn <- tryCatch(ProtConnCLa(nodes = nodes,
                                       region = Ecoreg_sel,
                                       thintersect = thintersect,
                                       attribute = attribute,
                                       distance = distance,
                                       transboundary = transboundary,
                                       distance_thresholds = distance_thresholds,
                                       probability = probability,
                                       plot = FALSE, intern = FALSE),  error = function(err)err)

      if (inherits(protconn, "error")){
        stop(protconn)
      }

      Ecoreglist <- tryCatch(lapply(protconn, function(x){
        n <- as.vector(x[[1]][[3]])
        x2 <- t(x[[1]][[4]]) %>% as.data.frame()
        colnames(x2) <- n
        return(x2)}), error = function(err)err)

      if (inherits( Ecoreglist, "error")){
        Ecoreglist <- lapply(protconn, function(x){
          n <- as.vector(x[[1]][[3]])
          x2 <- t(x[[1]][[4]]) %>% as.data.frame()
          colnames(x2) <- n
          return(x2)})
      }

      ECAlist <- tryCatch(lapply(protconn, function(x){
        n <- as.vector(x[[1]][[1]])
        x2 <- t(x[[1]][[2]]) %>% as.data.frame()
        colnames(x2) <- n
        return(x2)}), error = function(err)err)

      if (inherits(ECAlist, "error")){
        ECAlist <- lapply(protconn, function(x){
          n <- as.vector(x[[1]][[1]])
          x2 <- t(x[[1]][[2]]) %>% as.data.frame()
          colnames(x2) <- n
          return(x2)})
      }
      result <- list(Ecoreglist, ECAlist)
      return(result)}, .progress = intern), error = function(err)err)
    }

  if (inherits(protconn_result, "error")){
      setwd(ttt.2)
      stop(protconn_result)
      } else {
        Ecoreglist <- lapply(protconn_result, function(x){x[[1]]})
        ECAlist <- lapply(protconn_result, function(x){x[[2]]})
        data_1 <- list()
        for(i in 1:length(distance_thresholds)){
          data_2 <- list()
          for(j in 1:length(regions$ID_Temp)){
            data_2[[j]] <- Ecoreglist[[j]][[i]]
            }
          ProtConnEcor <- do.call(rbind, data_2)
          rownames(ProtConnEcor) <- NULL
          ProtConnEcor$ID_Temp <- regions$ID_Temp
          ####ECA
          ECA_2 <- list()
          for(j in 1:length(regions$ID_Temp)){
            ECA_2[[j]] <- ECAlist[[j]][[i]]
            }
      ECA_3 <- do.call(rbind, ECA_2)
      rownames(ECA_3) <- NULL
      ECA_3 <- ECA_3[,c(1:2)]
      ECA_3$ID_Temp <- regions$ID_Temp
      ####
      ProtConnEcor <- base::merge(ECA_3, ProtConnEcor, by = "ID_Temp")
      regions.2 <- base::merge(regions, ProtConnEcor, by = "ID_Temp")
      regions.2$ID_Temp <- NULL
      results <- list()
      results[[2]] <- regions.2
      DataProtconn <- regions.2
      sf::st_geometry(DataProtconn) <- NULL
      #
      DataProtconn_2 <-  tryCatch(ProtConnStat(DataProtconn, ci = CI, nr = 1000), error = function(err)err)
      if(inherits(DataProtconn_2, "error")){
        DataProtconn_2 <-  tryCatch(ProtConnStat(DataProtconn, ci = NULL), error = function(err)err)
        if(inherits(DataProtconn_2, "error")){
          DataProtconn_2 <- DataProtconn
        }
      }
      #
      DataProtconn_2$Indicator <- rownames(DataProtconn_2)
      DataProtconn_2 <- DataProtconn_2[moveme(names(DataProtconn_2), "Indicator first")]
      rownames(DataProtconn_2) <- NULL
      names(DataProtconn_2)[1:2] <- c("ProtConn indicator", "Values(%)")
      DataProtconn_2[c(2:ncol(DataProtconn_2))] <- round(DataProtconn_2[c(2:ncol(DataProtconn_2))], 3)

      if(ncol(DataProtconn_2) > 2){
        DataProtconn <- formattable::formattable(DataProtconn_2, align = c("l", rep("r", NCOL(DataProtconn_2) - 1)),
                                    list(`ProtConn indicator` = formattable::formatter("span", style = ~ formattable::style(color = "#636363", font.weight = "bold")),
                                         `Values(%)` = formattable::color_tile("white", "#F88B13"),
                                         formattable::area(col = 3:4) ~ formattable::color_tile("white", "#CE5D9B")))
      } else {
        DataProtconn <- formattable::formattable(DataProtconn_2, align = c("l", rep("r", NCOL(DataProtconn_2) - 1)),
                                    list(`ProtConn indicator` = formattable::formatter("span", style = ~ formattable::style(color = "#636363", font.weight = "bold")),
                                         `Values(%)` = formattable::color_tile("white", "#F88B13")))
      }

      results[[1]] <- DataProtconn
      names(results) <- c(paste0("ProtConn_", distance_thresholds[i]),
                          paste0("ProtConn_overall_", distance_thresholds[i]))
      ####################
      if(isTRUE(plot)){
        plots <- list()
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

          plot_protconn1 <- ggplot2::ggplot(dacc, ggplot2::aes(x = name, y = Values, fill = name)) +
            ggplot2::geom_bar(position = ggplot2::position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
            ggplot2::geom_errorbar(position = ggplot2::position_dodge(), width = 0.2, ggplot2::aes(ymin = min, ymax = max)) +
            ggplot2::labs(title = paste0("ProtConn Indicators: ", distance_thresholds[i]), x = "", y = "Percentage (%)", size = ggplot2::rel(1.2)) +
            ggplot2::theme_bw()  +
            ggplot2::theme(plot.title = ggplot2::element_text(color = "#252525", size = ggplot2::rel(1.4), hjust = 0.5, face = "bold"),
                  axis.title= ggplot2::element_text(color = "#252525", size = ggplot2::rel(1.2)),
                  legend.title= ggplot2::element_blank(),
                  legend.text = ggplot2::element_text(colour = "#252525", size = ggplot2::rel(1.2)),
                  axis.text= ggplot2::element_text(colour = "#525252", size = ggplot2::rel(1)))+
            ggplot2::scale_fill_manual(values = dacc$col) +
            ggplot2::geom_hline(ggplot2::aes(yintercept = 17, linetype = "Aichi Target (17%)"), colour = 'black', size = 1.2) +
            ggplot2::scale_linetype_manual(name = " Aichi Target", values = c(2, 2),
                                  guide = ggplot2::guide_legend(override.aes = list(color = c("black"), size = 0.8)))
          plots[[1]] <- plot_protconn1
        }

        dacc <- as.data.frame(DataProtconn[c(7:10), c(1:3)])
        dacc$name <- c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]")
        dacc$Indicator <- NULL
        dacc$name <- factor(dacc$name, levels = c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]"))
        dacc$col <- c("#253494", "#2c7fb8", "#41b6c4", "#7fcdbb")
        names(dacc)[2] <- "Values"
        dacc <- dacc[which(dacc$Values > 0), ]
        if(nrow(dacc) > 1){
          dacc$max <- dacc$Values+dacc$SD
          dacc$min <- dacc$Values-dacc$SD
          dacc[6:7] <- apply(dacc[6:7], 2, function(x){
            x[which(x >100)]<- 100
            x[which(x <0)]<- 0
            return(x)})
          plot_protconn2 <- ggplot2::ggplot(dacc, aes(x = name, y = Values, fill = name)) +
            ggplot2::geom_bar(position = ggplot2::position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
            ggplot2::geom_errorbar(position = ggplot2::position_dodge(), width = 0.2, ggplot2::aes(ymin = min, ymax = max)) +
            ggplot2::labs(title = "Protected connected fraction", x = "", y = "Percentage (%)", size = ggplot2::rel(1.2)) +
            ggplot2::theme_bw()  +
            ggplot2::theme(plot.title = ggplot2::element_text(color = "#252525", size = ggplot2::rel(1.4), hjust = 0.45, face = "bold"),
                  axis.title= ggplot2::element_text(color = "#252525", size = ggplot2::rel(1.2)),
                  plot.margin = ggplot2::margin(0, 5.3, 0, 0.3, "cm"),
                  legend.title= ggplot2::element_blank(),
                  legend.text = ggplot2::element_text(colour = "#252525", size = ggplot2::rel(1.2)),
                  axis.text= ggplot2::element_text(colour = "#525252", size = ggplot2::rel(1)))+
            ggplot2::scale_fill_manual(values = dacc$col)
          plots[[2]] <- plot_protconn2
        }
        plots <- plyr::compact(plots)

        if(length(plots) == 2){
          figure <- ggpubr::ggarrange(plots[[1]], plots[[2]],
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
        sf::write_sf(regions.2, paste0(write, "ProtConn_", distance_thresholds[i], ".shp"), delete_layer = T)
        write.csv(DataProtconn, paste0(write, "SummaryStats_", distance_thresholds[i], ".csv"))
        if(isTRUE(plot)){
          if(!is.character(figure)){
            tiff(paste0(write, "ProtConn_plot_d", distance_thresholds[i], ".tif"), width = 806, height = 641)
            print(figure)
            dev.off()
          }
        }
      }
      data_1[[i]] <- results
    }
        setwd(ttt.2)
        names(data_1) <- paste0("ProtConn_", distance_thresholds)
        return(data_1)
    }
  }
