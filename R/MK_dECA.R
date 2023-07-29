#' Estimate the ECA, rECA, dA and dECA.
#'
#' Equivalent Connected Area (ECA; if the area is used as attribute) or Equivalent Connectivity index (EC)
#' @param nodes list of objects class sf, SpatialPolygonsDataFrame or raster. Nodes of each time to analyze.
#' The shapefiles must be in a projected coordinate system.
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute in the data table selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit, "Makurhini::unit_covert()"
#' compatible unit(e.g., "m2", "km2", "ha"). Default equal to hectares "m2".
#' @param distance list. Distance parameters. For example: type, resistance,or keep. For "type" choose one of the
#'  distances: "centroid" (faster), "edge", "least-cost" or "commute-time". If the type is equal to "least-cost"
#'  or "commute-time", then you have to use the "resistance" argument. You can place a list with resistances
#'  where there must be one resistance for each time/scenario of patches in the **nodes** parameter, for example,
#'  if nodes has a list with two time patches then you can use two resistances one for time 1 and one for time 2:
#'  distance(type = "least-cost", resistance = list(resistanceT1, resistanceT2)). To See more arguments consult
#'  the help function of distancefile().
#' @param metric character. Choose a connectivity metric: "IIC" considering topologycal distances or "PC"
#' considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g.,
#' 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' If probability = NULL, then it will be the inverse of the mean dispersal distance for the species (1/α; Hanski and Ovaskainen 2000).
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters).
#' For example, distance_threshold = 30000; two or more specific distances:
#' distance_threshold = c(30000, 50000); sequence distances: distance_threshold = seq(10000,100000, 10000).
#' @param LA numeric. Maximum landscape attribute ("units" equal to "area_unit", default equal to "ha").
#' @param plot logical. Also, you can provide the corresponding year for each period of
#' time analyzed, e.g., c("2011", "2014", "2017")
#' @param parallel numeric. Specify the number of cores to use for parallel processing, default = NULL.
#' Parallelize the function using furrr package and multiprocess plan.
#' @param write character. Path and name of the output ".csv" file
#' @param intern logical. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @return Table with:\cr\cr
#' -  Time: name of the time periods, name of the model or scenario (are taken from the name of the elements of the list of nodes or the plot argument)\cr
#' -  Max. Landscape atrtribute: maximum landscape attribute\cr
#' -  Habitat area,\cr
#' -  Distance threshold: it is usually a dispersal threshold associated with one or many species and it is set by the user.\cr
#' -  ECA: Equivalent Connected Area or Equivalent Connectivity
#' -  Normalized_ECA (% of LA): relative connectivity (percentage)\cr
#' -  Normalized_ECA (% of habitat area): relative connectivity (percentage)\cr
#' -  dA: delta Area between times (percentage)\cr
#' -  dECA: delta ECA between times (percentage)\cr
#' -  rECA: relativized ECA (dECA/dA). According to Liang et al. (2021) "an rECA value greater than 1 indicates that habitat changes result in a
#' disproportionately large change in habitat connectivity, while a value lower than 1 indicates connectivity
#' changes due to random habitat changes (Saura et al. 2011; Dilts et al. 2016)".\cr
#' -  dA/dECA comparisons: comparisons between dA and dECA\cr
#' -  Type of change: Type of change using the dECAfun() and the difference between dA and dECA.\cr
#' @references \url{www.conefor.org}\cr\cr
#' Saura, S., Estreguil, C., Mouton, C., & Rodríguez-Freire, M. (2011). Network analysis to assess landscape connectivity trends: Application to European forests (1990-2000). Ecological Indicators, 11(2), 407–416.
#' https://doi.org/10.1016/j.ecolind.2010.06.011 \cr\cr
#' Herrera, L. P., Sabatino, M. C., Jaimes, F. R., & Saura, S. (2017). Landscape connectivity and the role of small habitat patches as stepping stones: an assessment of the grassland biome in South America. Biodiversity and Conservation, 26(14), 3465–3479.
#' https://doi.org/10.1007/s10531-017-1416-7\cr\cr
#' Liang, J., Ding, Z., Jiang, Z., Yang, X., Xiao, R., Singh, P. B., ... & Hu, H. (2021). Climate change, habitat connectivity, and conservation gaps: a case study of four ungulate species endemic to the Tibetan Plateau. Landscape Ecology, 36(4), 1071-1087.
#' Dilts TE, Weisberg PJ, Leitner P, Matocq MD, Inman RD, Nussear KE, Esque TC (2016) Multi-scale connectivity and graph theory highlight critical areas for conservation under climate change. Ecol Appl 26:1223–1237
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(sf)
#'
#' data("list_forest_patches", package = "Makurhini")
#' data("study_area", package = "Makurhini")
#' class(list_forest_patches)
#'
#' Max_attribute <- unit_convert(st_area(study_area), "m2", "ha")
#'
#' dECA_test <- MK_dECA(nodes= list_forest_patches, attribute = NULL, area_unit = "ha",
#'                   distance = list(type= "centroid"), metric = "PC",
#'                   probability = 0.05, distance_thresholds = 5000,
#'                   LA = Max_attribute, plot= c("1993", "2003", "2007", "2011"),
#'                   intern = TRUE)
#' dECA_test
#'
#' }
#' @export
#' @importFrom magrittr %>%
#' @importFrom purrr compact map_dfr
#' @importFrom methods as
#' @importFrom formattable formattable color_bar proportion formatter percent style
#' @importFrom ggplot2 ggplot geom_bar aes geom_text theme element_blank element_text labs ggtitle scale_fill_manual element_line ggsave
#' @importFrom utils write.csv txtProgressBar setTxtProgressBar installed.packages
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map

MK_dECA <- function(nodes,
                    attribute = NULL,
                    area_unit = "m2",
                    distance = list(type = "centroid", resistance = NULL),
                    metric = "IIC",
                    probability = NULL,
                    distance_thresholds = NULL,
                    LA = NULL,
                    plot = FALSE, parallel = NULL,
                    write = NULL, intern = TRUE){
  . = NULL
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

  if (is.null(distance_thresholds)) {
    stop("error missing numeric distance threshold(s)")
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }

    if(!is.character(write)){
      write = NULL
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

  options(warn = -1); listT <- compact(nodes)

  if(distance$type == "least-cost" | distance$type == "commute-time"){
    if(is.null(distance$resistance)){
      stop(paste0("error missing resistance raster, e.g., resistance(type =", distance$type,
                  ", resistance = MISSING_RESISTANCE)"))
    } else {
      if(class(distance$resistance)[1] == "list"){
        if(length(distance$resistance) > 1){
          if(length(distance$resistance) != length(listT)){
            stop(paste0("the length of the node list (", length(listT), ") must be equal to the length of resistance list (", length(distance$resistance), ")"))
          }
        }
      }
    }
  } else {
    dist_param <- distance
  }

  if(class(listT[[1]])[1] != "RasterLayer"){
    listT <- lapply(listT, function(x) { if(class(x)[1] != "sf") {
      x <- st_as_sf(x); x$IdTemp <- 1:nrow(x)
    } else {
      x$IdTemp <- 1:nrow(x)
    }
      return(x)
    })
  }
  #
  time <- as.vector(1:length(listT)) %>% as.character()
  #
  if(class(listT[[1]])[1] != "RasterLayer"){
    DECA <- map_dfr(listT, function(x){
      x.1 <- unit_convert(sum(st_area(x, byid = T)), "m2", area_unit) %>%
        as.data.frame()
      return(x.1)
      }); id = "IdTemp"
  } else {
    nres <- unit_convert(res(listT[[1]])[1]^2, "m2", area_unit)
    DECA <- lapply(listT, function(x){
      x1 <- as.data.frame(table(x[])); x2 <- sum(x1$Freq * nres)
      return(x2)}) %>% do.call(rbind, .) %>% as.data.frame(as.numeric(.)); id = NULL
  }

  DECA <- cbind(time, LA, DECA); rownames(DECA) <- NULL; colnames(DECA)[3]<-"Area"

  #ECA
  x = NULL; y = NULL

  if(isTRUE(intern)){
    pb <- txtProgressBar(0, length(listT), style = 3)
  }

  if(is.null(parallel)){
    ECA <- tryCatch(lapply(1:length(listT), function(x){
      x.1 <- listT[[x]]

      if(isTRUE(intern)){
        setTxtProgressBar(pb, x)
      }

      if(nrow(x.1) < 2){
        x.1 <- unit_convert(sum(st_area(x.1, byid = T)), "m2", area_unit)
        ECA_metric <- map_dfr(distance_thresholds, function(y) {
          tab1 <- if((100 * (x.1 / LA)) >= 100){LA}else{x.1}; tab1 <- data.frame(Value = tab1)
          return(tab1)})
        ECA_metric$Distance <- distance_thresholds; names(ECA_metric)[1] <- "ECA"
      } else {
        ECA_metric <- map_dfr(distance_thresholds, function(y) {
          if(distance$type == "least-cost" | distance$type == "commute-time"){
            if(class(distance$resistance)[1] == "list"){
              dist_param <- distance; dist_param$resistance <- distance$resistance[[x]]
            } else {
              dist_param <- distance
            }
          } else {
            dist_param <- distance
          }

          tab1 <- MK_dPCIIC(nodes = x.1, attribute = attribute,
                            restoration = NULL,
                            distance = dist_param, area_unit = area_unit,
                            metric = metric, probability = probability,
                            distance_thresholds = y,
                            overall = TRUE, onlyoverall = TRUE,
                            LA = LA, rasterparallel = FALSE, write = NULL)
          return(tab1[2,])
        })
        ECA_metric$Distance <- distance_thresholds; names(ECA_metric)[2] <- "ECA"; ECA_metric <- ECA_metric[,2:3]
      }
      return(ECA_metric)
    }), error = function(err) err)
  } else {
    works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}
    if(.Platform$OS.type == "unix") {
      strat <- future::multicore
    } else {
      strat <- future::multisession
    }
    plan(strategy = strat, gc = TRUE, workers = works)
    ECA <- tryCatch(future_map(1:length(listT), function(x) {
      x.1 <- listT[[x]]

      if(nrow(x.1) < 2){
        x.1 <- unit_convert(sum(st_area(x.1, byid = T)), "m2", area_unit)
        ECA_metric <- map_dfr(distance_thresholds, function(y) {
          tab1 <- if((100 * (x.1 / LA)) >= 100){LA}else{x.1}; tab1 <- data.frame(Value = tab1)
          return(tab1)})
        ECA_metric$Distance <- distance_thresholds; names(ECA_metric)[1] <- "ECA"
      } else {
        ECA_metric <- map_dfr(distance_thresholds, function(y) {
          if(distance$type == "least-cost" | distance$type == "commute-time"){
            if(class(distance$resistance)[1] == "list"){
              dist_param <- distance; dist_param$resistance <- distance$resistance[[x]]
            } else {
              dist_param <- distance
            }
          } else {
            dist_param <- distance
          }

          tab1 <- MK_dPCIIC(nodes = x.1, attribute = attribute,
                            restoration = NULL,
                            distance = dist_param, area_unit = area_unit,
                            metric = metric, probability = probability,
                            distance_thresholds = y,
                            overall = TRUE, onlyoverall = TRUE,
                            LA = LA, rasterparallel = FALSE, write = NULL)
          return(tab1[2,])
        })
        ECA_metric$Distance <- distance_thresholds; names(ECA_metric)[2] <- "ECA"; ECA_metric <- ECA_metric[,2:3]
      }
      return(ECA_metric)
    }, .progress = TRUE),  error = function(err) err)
    close_multiprocess(works)
  }

  if(inherits(ECA, "error")){
    stop("review ECA parameters: nodes, distance file or LA")
  } else {
    ECA2 <- lapply(distance_thresholds, function(x){
      Tab_ECA <- map_dfr(ECA, function(y){ y[which(y$Distance == x),] })
      DECA.2 <- cbind(DECA, Tab_ECA)

      AO <- LA; dArea <- (DECA.2$Area *100)/LA
      DECA.2$Normalized_ECA1 <- (DECA.2$ECA* dArea)/DECA.2$Area
      DECA.2$Normalized_ECA2 <- (DECA.2$ECA*100)/DECA.2$Area
      AO.2 <- DECA.2$Area; ECA.2 <- cbind(DECA.2[, 4])
      DECA.2$dA[1] <- (((DECA.2$Area[1] - AO)/AO) * 100)
      DECA.2$dA[2:nrow(DECA.2)] <- ((DECA.2$Area[2:nrow(DECA.2)] - AO.2)/AO.2) *100

      DECA.2$dECA[1] <- (((DECA.2$ECA[1] - AO)/AO) * 100)
      DECA.2$dECA[2:nrow(DECA.2)] <- ((DECA.2$ECA[2:nrow(DECA.2)] - ECA.2)/ECA.2) * 100

      DECA.2[,2:ncol(DECA.2)] <- round(DECA.2[,2:ncol(DECA.2)], 3)

      DECA.3 <- data.frame(time = DECA.2$time, Type_Change = "", Type ="")

      for(i in 1:nrow(DECA.2)){
        DECA.3[i,2] <- dECAfun(DECA.2$dECA[i], DECA.2$dA[i])
        DECA.3[i,3] <- dECAfun2(DECA.2$dECA[i], DECA.2$dA[i])
      }

      names(DECA.3)[2:3] <- c("dA/dECA comparisons", "Type of change")
      DECA.2 <- DECA.2[,c(1:3, 5, 4, 6:ncol(DECA.2))]
      DECA.2$rECA <- DECA.2$dECA/DECA.2$dA
      DECA.4 <- cbind(DECA.2, DECA.3); DECA.4[11] <- NULL

      names(DECA.4)[c(1:7)] <- c("Time",
                                 paste0("Max. Landscape attribute (",area_unit,")"),
                                 paste0("Habitat area (",area_unit,")"),
                                 "Distance threshold",
                                 paste0("ECA (",area_unit,")"),
                                 "Normalized ECA (% of LA)",
                                 "Normalized ECA (% of habitat area)"                                 )

      if(is.character(plot) & length(plot) == nrow(DECA.4)){
        DECA.4$Time <- factor(plot, levels = plot)
      } else {
        if(is.null(names(listT))){
          DECA.4$Time <- rownames(DECA.4)
        } else {
          DECA.4$Time <- names(listT)
        }
      }

      rownames(DECA.4) <- NULL

      DECA.4 <- formattable(DECA.4, align = c("l", rep("r", NCOL(DECA.4) - 1)),
                            list(`Habitat area (ha)`= color_bar("#94D8B1", proportion),
                                 `ECA (ha)` = color_tile("#FAE9EA", "#E8878A"),
                                 `Normalized ECA (% of LA)` = color_tile("#FFEDCC", "orange"),
                                 `Normalized ECA (% of habitat area)` = color_tile("#FFEDCC", "orange"),
                                 `dA` = formatter("span",style = ~ style(color = ifelse(`dA` > 0, "forestgreen", "red"))),
                                 `dECA` = formatter("span",style = ~ style(color = ifelse(`dECA` > 0, "forestgreen", "red"))),
                                 `rECA` = formatter("span",style = ~ style(color = ifelse(`rECA` > 0, "#404040", "red")))))
      return(DECA.4)})
    #
    if (!is.null(write)){
      write.csv(do.call(rbind, ECA2), write, row.names = FALSE)
    }

    if(isTRUE("ggplot2" %in% rownames(installed.packages()))){
      if(isTRUE(plot) | is.character(plot)){
        if(isTRUE(plot)){
          plot <- paste0("Time", 1:length(nodes))
        }
        plot <- factor(plot, levels = plot)

        ECAplot <- lapply(ECA2, function(x){
          ECA4 <- (x[[3]] * 100)/ LA; Loss <- 100 - ECA4; ECA4 <- cbind(ECA4, Loss) %>% as.data.frame()
          names(ECA4)[1] <- "Habitat"; ECA4$"Connected habitat" <- x[[7]]; ECA4$Year <- plot

          #Table 1
          ECA5 <- ECA4
          ECA4 <- data.frame(Year = rep(ECA4$Year,3),
                             variable = c(rep("Habitat", length(ECA4$Year)),
                                          rep("Loss", length(ECA4$Year)),
                                          rep("Connected habitat", length(ECA4$Year))),
                             percentage = c(ECA4$Habitat, ECA4$Loss, ECA4$`Connected habitat`),
                             check.names = FALSE)

          ECA4$variable <- factor(ECA4$variable, levels = c("Loss", "Habitat", "Connected habitat"))

          #Table 2
          ECA5$`Connected habitat` <- (ECA5$`Connected habitat` * ECA5$Habitat)/100
          ECA5$Habitat <- ECA5$Habitat - ECA5$`Connected habitat`

          ECA5 <- data.frame(Year = rep(ECA5$Year,3),
                             variable = c(rep("Habitat", length(ECA5$Year)),
                                          rep("Loss", length(ECA5$Year)),
                                          rep("Connected habitat", length(ECA5$Year))),
                             percentage = c(ECA5$Habitat, ECA5$Loss, ECA5$`Connected habitat`),
                             check.names = FALSE)
          ECA5$variable <- factor(ECA5$variable, levels = c("Loss", "Habitat", "Connected habitat"))
          ECA5$Text <- ECA4$percentage
          ECA5$Text[which(ECA5$variable == "Connected habitat")] <- ECA5$percentage[which(ECA5$variable == "Connected habitat")]

          ECA5$pos <- ECA5$percentage/2
          ECA5$pos[which(ECA5$variable == "Habitat")] <- (ECA5$percentage[which(ECA5$variable == "Habitat")]/2) + ECA5$percentage[which(ECA5$variable == "Connected habitat")]
          ECA5$pos[which(ECA5$variable == "Loss")] <- (ECA5$percentage[which(ECA5$variable == "Loss")]/2) + ECA5$Text[which(ECA5$variable == "Habitat")]
          #Plot
          pcolors <- c("#443A82", "#30678D", "#26AE7F")

          if (distance$type  %in% c("centroid", "edge")){
            if(is.null(distance$distance_unit)){
              dp_text <-  paste0(unique(x$Distance), " m")
            } else {
              dp_text <-  paste0(unique(x$Distance)," " ,distance$distance_unit)
            }
          } else {
            dp_text <-  paste0(unique(x$Distance), " cost-weighted distance")
          }

          p4 <- ggplot() +
            geom_bar(aes(y = ECA5$percentage, x = ECA5$Year, fill = ECA5$variable), data = ECA5, stat="identity") +
            geom_text(data = ECA5, aes(x = ECA5$Year, y = ECA5$pos, label = paste0(round(ECA5$Text, 2), "%")),
                      colour = "white", family = "Tahoma", size = 4, fontface = "bold") +
            theme(legend.position = "bottom", legend.direction = "horizontal",
                  legend.title = element_blank(),
                  legend.text = element_text(size = 11)) +
            labs(x = "Time", y = "% of landscape") +
            ggtitle(paste0("ECA: Dispersal distance = ", dp_text)) +
            scale_fill_manual(values = pcolors) +
            theme(axis.line = element_line(size = 1, colour = "black"),
                  panel.grid.major = element_line(colour = "#d3d3d3"),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank()) +
            theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
                  text = element_text(family = "Tahoma", size = 12),
                  axis.text.x = element_text(colour = "black", size = 10),
                  axis.text.y = element_text(colour = "black", size = 10))

          if(!is.null(write)){
            ggsave(paste0(write, "_", unique(x$Distance), '_ECA.tif'), plot = p4, device = "tiff", width = 14,
                   height = 10, compression = "lzw", dpi = "retina", scale = 0.7)}
          return(p4)})
        ECA_result <- list()
        for (i in 1:length(ECA2)){
          ECA_result[[i]] <- list(ECA2[[i]], ECAplot[[i]])
        }

        names(ECA_result) <- paste(distance_thresholds)
        ECA2 <- ECA_result
      }
    } else {
      message("You need install ggplot2 to plot dECA")
      plot = FALSE
    }


    if(length(distance_thresholds) == 1){
      if((isTRUE(plot) | is.character(plot))){
        print(ECA2[[1]][[2]]); ECA2 <- ECA2[[1]][[1]]
      } else {
        ECA2 <- ECA2[[1]]
      }
    }
  }
  return(ECA2)
}
