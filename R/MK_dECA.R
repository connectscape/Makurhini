#' ECA, dA and dECA.
#'
#' Equivalent Connected Area (ECA; if the area is used as attribute) or Equivalent Connectivity index (EC)
#' @param nodes list of objects of class sf, SpatialPolygonsDataFrame or raster. Nodes of each period of time to analyze.
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute in the data table selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit, "Makurhini::unit_covert()"
#' compatible unit(e.g., "m2", "km2", "ha"). Default equal to hectares "ha".
#' @param distance list. Distance parameters. For example: type, resistance,or keep. For "type" choose one of the
#'  distances: "centroid" (faster), "edge", "least-cost" or "commute-time". If the type is equal to "least-cost"
#'  or "commute-time", then you have to use the "resistance" argument. To See more arguments consult
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
#' @param parallel logical. Parallelize the function using furrr package and multiprocess plan, default = FALSE.
#' @param write character. Path and name of the output ".csv" file
#' @return Table with:\cr
#' A: Area in km2\cr
#' ECA: ECA value\cr
#' Normalized_ECA: Relative connectivity (percentage)\cr
#' dA: Delta Area between times (percentage)\cr
#' dECA: Delta ECA between times (percentage)\cr
#' Type_change: Type of change using the dECAfun() and the difference between dA and dECA.\cr
#' @references \url{www.conefor.org}\cr
#' Saura, S., Estreguil, C., Mouton, C., & Rodríguez-Freire, M. (2011). Network analysis to assess landscape connectivity trends: Application to European forests (1990-2000). Ecological Indicators, 11(2), 407–416.
#' https://doi.org/10.1016/j.ecolind.2010.06.011 \cr
#' Herrera, L. P., Sabatino, M. C., Jaimes, F. R., & Saura, S. (2017). Landscape connectivity and the role of small habitat patches as stepping stones: an assessment of the grassland biome in South America. Biodiversity and Conservation, 26(14), 3465–3479.
#' https://doi.org/10.1007/s10531-017-1416-7
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(rgeos)
#'
#' data("list_forest_patches", package = "Makurhini")
#' data("study_area", package = "Makurhini")
#' class(list_forest_patches)
#'
#' Max_attribute <- unit_convert(gArea(study_area), "m2", "ha")
#'
#' dECA_test <- MK_dECA(nodes= list_forest_patches, attribute = NULL, area_unit = "ha",
#'                   distance = list(type= "centroid"), metric = "PC",
#'                   probability = 0.05, distance_thresholds = 5000,
#'                   LA = Max_attribute, plot= c("1993", "2003", "2007", "2011"))
#' dECA_test
#'
#' }
#' @export
#' @importFrom magrittr %>%
#' @importFrom purrr compact map map_dfr
#' @importFrom methods as
#' @importFrom rgeos gArea
#' @importFrom dplyr progress_estimated summarize
#' @importFrom plyr ddply .
#' @importFrom formattable formattable color_bar proportion formatter percent style
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_bar aes geom_text theme element_blank element_text labs ggtitle scale_fill_manual element_line ggsave
#' @importFrom utils write.csv
#' @importFrom future multiprocess plan availableCores
#' @importFrom furrr future_map future_map_dfr
#' @importFrom rlang .data

MK_dECA <- function(nodes,
                    attribute = NULL,
                    area_unit = "ha",
                    distance = list(type = "centroid", resistance = NULL),
                    metric = "IIC",
                    probability = NULL,
                    distance_thresholds = NULL,
                    LA = NULL,
                    plot = FALSE, parallel = FALSE,
                    write = NULL){
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

  options(warn = -1)
  listT <- compact(nodes)

  if(class(listT[[1]])[1] != "RasterLayer"){
    listT <- map(listT, function(x) { if(class(x)[1] == "sf") {
      x <- as(x, 'Spatial')
      x@data$IdTemp <- 1:nrow(x)
    } else {
      x@data$IdTemp <- 1:nrow(x)
    }
      return(x)
    })
  }

  #
  scenary <- as.vector(1:length(listT)) %>% as.character()
  #
  if(class(listT[[1]])[1] == "SpatialPolygonsDataFrame"){
    DECA <- map_dfr(listT, function(x){unit_convert(sum(gArea(x, byid = T)), "m2", area_unit) %>% as.data.frame()})
    id = "IdTemp"
  } else {
    nres <- unit_convert(res(listT[[1]])[1]^2, "m2", area_unit)
    DECA <- map(listT, function(x){
      x1 <- as.data.frame(table(x[]))
      sum(x1$Freq * nres)})
    DECA <- do.call(rbind, DECA) %>% as.data.frame(as.numeric(.))
    id = NULL
  }

  #DECA <- do.call(rbind, DECA) %>% as.data.frame(as.numeric(.))
  DECA <- cbind(scenary, DECA)
  rownames(DECA) <- NULL
  colnames(DECA)[2]<-"Area"

  #ECA
  x = NULL
  y = NULL
  if(isFALSE(parallel)){
    pb <- progress_estimated(length(listT), 0)
    ECA <- tryCatch(map(listT, function(x){
      pb$tick()$print()

      ECA_metric <-  map_dfr(distance_thresholds, function(y) {
        tab1 <- MK_dPCIIC(nodes = x, attribute = attribute,
                          restoration = NULL,
                          distance = distance, area_unit = area_unit,
                          metric = metric, probability = probability,
                          distance_thresholds = y,
                          overall = TRUE, onlyoverall = TRUE,
                          LA = LA, rasterparallel = FALSE, write = NULL)
        tab1 <- tab1[2,2]
        return(tab1)
      })

      #ECA_metric2 <- do.call(rbind,  ECA_metric)
      ECA_metric2 <- cbind(ECA_metric, distance_thresholds)
      ECA_metric2 <- as.data.frame(ECA_metric2)
      names(ECA_metric2) <- c("ECA", "Distance")
      return(ECA_metric2)
    }), error = function(err) err)
  } else {
    works = as.numeric(availableCores())-
      plan(strategy = multiprocess, gc = TRUE, workers = works)
    ECA <- tryCatch(future_map(listT, function(x) {
      ECA_metric <-  future_map_dfr(distance_thresholds, function(y) {
        tab1 <- MK_dPCIIC(nodes = x, attribute = attribute,
                          restoration = NULL,
                          distance = distance, area_unit = area_unit,
                          metric = metric, probability = probability,
                          distance_thresholds = y,
                          overall = TRUE, onlyoverall = TRUE,
                          LA = LA, rasterparallel = FALSE, write = NULL)
        tab1 <- tab1[2,2]
        return(tab1)
      })
      ECA_metric2 <- cbind(ECA_metric, distance_thresholds)
      ECA_metric2 <- as.data.frame(ECA_metric2)
      names(ECA_metric2) <- c("ECA", "Distance")
      return(ECA_metric2)
    }, .progress = TRUE),  error = function(err) err)
    close_multiprocess(works)
  }

  if(inherits(ECA, "error")){
    stop("review ECA parameters: nodes, distance file or LA")
  } else {

    ECA2 <- map(distance_thresholds, function(x){
      Tab_ECA <- map(ECA, function(y){ y[which(y$Distance == x),] })
      Tab_ECA <- do.call(rbind, Tab_ECA)
      Tab_ECA <- cbind(DECA, Tab_ECA)
      return(Tab_ECA)})

    ECA3 <- map(ECA2, function(x){
      DECA.2 <- x
      AO <- LA

      DECA.2$Normalized_ECA <- (DECA.2$ECA*100)/DECA.2$Area
      AO.2 <- DECA.2$Area
      ECA.2 <- cbind(DECA.2[1:nrow(DECA.2), 3])

      DECA.2$dA[1] <- (((DECA.2$Area[1] - AO)/AO) * 100)
      DECA.2$dA[2:nrow(DECA.2)] <- ((DECA.2$Area[2:nrow(DECA.2)] - AO.2)/AO.2) *100

      DECA.2$dECA[1] <- (((DECA.2$ECA[1] - AO)/AO) * 100)
      DECA.2$dECA[2:nrow(DECA.2)] <- ((DECA.2$ECA[2:nrow(DECA.2)] - ECA.2)/ECA.2) * 100

      DECA.2[,2:ncol(DECA.2)] <- round(DECA.2[,2:ncol(DECA.2)], 3)

      DECA.3 <- ddply(DECA.2, .(scenary), dplyr::summarize,
                      Type_Change = dECAfun(.data$dECA, .data$dA))
      DECA.3$Type <- ddply(DECA.2, .(scenary), dplyr::summarize,
                           Type = dECAfun2(.data$dECA, .data$dA))[[2]]
      names(DECA.3)[2:3] <- c("dA/dECA comparisons", "Type of change")
      DECA.4 <- cbind(DECA.2, DECA.3)
      DECA.4[8] <- NULL

      names(DECA.4)[c(1:3,5)] <- c("Scenary", paste0("Area (",area_unit,")"),
                                   paste0("ECA (",area_unit,")"), "Normalized ECA")

      if(is.character(plot) & length(plot) == nrow(DECA.4)){
        DECA.4$Scenary <- plot
      } else {
        DECA.4$Scenary <- rownames(DECA.4)
      }

      rownames(DECA.4) <- NULL

      DECA.4 <- formattable(DECA.4, align = c("l", rep("r", NCOL(DECA.4) - 1)),
                            list(`Area (ha)`= color_bar("#94D8B1", proportion),
                                 `Normalized ECA` = formatter("span", x ~ percent(x / 100)),
                                 `dA` = formatter("span",style = ~ style(color = ifelse(`dA` > 0, "green", "red"))),
                                 `dECA` = formatter("span",style = ~ style(color = ifelse(`dECA` > 0, "green", "red")))))
      return(DECA.4)
    })

    #
    if (!is.null(write)){
      write.csv(do.call(rbind, ECA3), write, row.names = FALSE)
    }
    ###plot
    if(isTRUE(plot) | is.character(plot)){
      if(isTRUE(plot)){
        plot = paste0("Time", 1:length(nodes))
      }

      ECAplot <- map(ECA3, function(x){
        ECA4 <- (x[2] * 100)/ LA
        names(ECA4) <- "Habitat"
        ECA4$Loss <- 100 - ECA4$Habitat
        ECA4$"Connected habitat" <- x[[5]]
        ECA4$Year <- plot

        #Table 1
        ECA5 <- ECA4
        ECA4 <- melt(ECA4, id="Year")
        names(ECA4)[3] <- "percentage"
        ECA4$variable <- factor(ECA4$variable, levels = c("Loss", "Habitat", "Connected habitat"))

        #Table 2
        ECA5$`Connected habitat` <- (ECA5$`Connected habitat` * ECA5$Habitat)/100
        ECA5$Habitat <- ECA5$Habitat - ECA5$`Connected habitat`
        ECA5 <- melt(ECA5, id="Year")
        names(ECA5)[3] <- "percentage"

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
          labs(x = "Year", y = "Percentage (%)") +
          ggtitle(paste0("Equivalent Connected Area: Dispersal distance = ", dp_text)) +
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
      for (i in 1:length(ECA3)){
        ECA_result[[i]] <- list(ECA3[[i]], ECAplot[[i]])
      }

      names(ECA_result) <- paste(distance_thresholds)
      ECA3 <- ECA_result
    }
    #
    if(length(distance_thresholds) == 1){
      if((isTRUE(plot) | is.character(plot))){
        ECA4 <- ECA3[[1]][[2]]
        ECA3 <- ECA3[[1]][[1]]
        print(ECA4)
      } else {
        ECA3 <- ECA3[[1]]
      }
    }
  }
  return(ECA3)
}
