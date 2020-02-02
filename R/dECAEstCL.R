#' ECA, dA and dECA estimation using the Conefor command line.
#'
#' Equivalent Connected Area (ECA; if the area is used as attribute) or Equivalent Connectivity index (EC)
#' @param nodes list of objects of class sf, sfc, sfg or SpatialPolygons. Nodes of each period of time to analyze.
#' @param attribute character. Column name with the attribute in the data table selected for the nodes.
#'  If NULL the node area will be used as node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit, udunits2 package compatible
#' unit (e.g., "km2", "cm2", "ha"). Default equal to hectares "ha".
#' @param distance list. Distance parameters. For example: type, resistance,or tolerance. For "type" choose one of the
#'  distances: "centroid" (faster), "edge", "least-cost" or "commute-time". If the type is equal to "least-cost"
#'  or "commute-time", then you have to use the "resistance" argument. To See more arguments consult
#'  the help function of distancefile().
#' @param metric character. Choose a connectivity metric: "IIC" considering topologycal distances or "PC"
#' considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g.,
#' 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections.
#' For example, distance_threshold = 30000; two or more specific distances:
#' distance_threshold = c(30000, 50000); sequence distances: distance_threshold = seq(10000,100000, 10000).
#' @param LA numeric. Maximum landscape attribute ("units" equal to "area_unit", default equal to "ha").
#' @param plot character. If it is not NULL, You have to provide the corresponding year for each period of
#' time analyzed, e.g., c("2011", "2014", "2017")
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
#' ruta <- system.file("extdata", "ECA_example.RData", package = "Makurhini")
#' load(ruta)
#' LA_area <- raster::area(LA)* 0.0001
#' dECA <- dECAEstCL(nodes= nodes, attribute = NULL, area_unit = "ha",
#'                   distance = list(type= "centroid"), metric = "IIC",
#'                   probability = NULL, distance_thresholds = 30000,
#'                   LA = LA_area, plot= c("1T", "2T", "3T", "4T"),
#'                   write = NULL)
#' dECA[[1]]
#' dECA[[2]]
#' }
#' @export
#' @importFrom magrittr %>%
#' @importFrom purrr compact map
#' @importFrom methods as
#' @importFrom udunits2 ud.convert
#' @importFrom rgeos gArea
#' @importFrom dplyr progress_estimated
#' @importFrom plyr ddply .
#' @importFrom dplyr summarize
#' @importFrom formattable formattable color_bar proportion formatter percent style
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom utils write.csv
dECAEstCL <- function(nodes,
                      attribute = NULL,
                      area_unit = "ha",
                      distance = list(type = "centroid", resistance = NULL),
                      metric = "IIC",
                      probability = NULL,
                      distance_thresholds = NULL,
                      LA = NULL,
                      plot = NULL,
                      write = NULL){
  options(warn = -1)
  listT <- compact(nodes)
  listT <- map(listT, function(x) { if(class(x)[1] == "sf") {
    x <- as(x, 'Spatial') } else { x }})

  #
  scenary <- as.vector(1:length(listT)) %>% as.character()
  DECA <- map(listT, function(x){ud.convert(sum(gArea(x, byid = T)), "m2", area_unit)})
  DECA <- do.call(rbind, DECA) %>% as.data.frame(as.numeric(.))
  DECA <- cbind(scenary, DECA)
  rownames(DECA) <- NULL
  colnames(DECA)[2]<-"Area"

  #ECA
  ttt.2 <- getwd()
  temp.1 <- paste0(tempdir(), "/TempInputs", sample(1:1000, 1, replace = T))
  dir.create(temp.1, recursive = T)
  setwd(temp.1)
  #############################
  pb <- progress_estimated(length(listT), 0)

  ECA <- tryCatch(map(listT, function(x){
    pb$tick()$print()
    x@data$IdTemp <- 1:nrow(x)
    nodesfile(x, id = "IdTemp", attribute, area_unit, write = paste0(temp.1,"/nodes.txt"))

    distancefile(x,  id = "IdTemp", type = distance$type, tolerance = distance$tolerance,
               resistance = distance$resistance, CostFun = distance$CostFun, ngh = distance$ngh,
               threshold = distance$threshold, mask = distance$mask,
               distance_unit = distance$distance_unit, distance$geometry_out,
               write = paste0(temp.1,"/Dist.txt"))

    if (is.null(distance$threshold)) {
      pairs = "all"
      } else {
        pairs = "notall"
        }

    ECA_metric <-  map(as.list(distance_thresholds), function(y){
      tab1 <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                         typeconnection = "dist", typepairs = pairs,
                         index = metric, thdist = y,
                         multdist = NULL, conprob = probability,
                         onlyoverall = TRUE, LA = LA,
                         nrestauration = FALSE,
                         prefix = NULL, write = NULL)
      tab1 <- tab1[[2]]
      tab1 <- tab1[[2,2]]
      return(tab1)})

    ECA_metric2 <- do.call(rbind,  ECA_metric)
    ECA_metric2 <- cbind(ECA_metric2, distance_thresholds)
    ECA_metric2 <- as.data.frame(ECA_metric2)
    names(ECA_metric2) <- c("ECA", "Distance")
    return(ECA_metric2)}), error = function(err)err)

  if(inherits(ECA, "error")){
    setwd(ttt.2)
    stop("review ECA parameters: nodes, distance file or LA")
  } else {
    ECA2 <- map(distance_thresholds, function(x){
      Tab_ECA <- map(ECA, function(y){ y[which(y$Distance == x),] })
      Tab_ECA <- do.call(rbind, Tab_ECA)
      Tab_ECA <- cbind(DECA, Tab_ECA)
      return(Tab_ECA)})
  #

  ECA3 <- map(ECA2, function(x){
      DECA.2 <- x
      AO <- LA
      DECA.2$Normalized_ECA <- (DECA.2$ECA*100)/DECA.2$Area
      AO.2 <- rbind(DECA.2[1:3, 2])
      ECA.2 <- cbind(DECA.2[1:3, 3])
      DECA.2$dA[1] <- (((DECA.2$Area[1] - AO)/AO) * 100)
      DECA.2$dA[2:4] <- ((DECA.2$Area[2:4] - AO.2)/AO.2) *100
      DECA.2$dECA[1] <- (((DECA.2$ECA[1] - AO)/AO) * 100)
      DECA.2$dECA[2:4] <- ((DECA.2$ECA[2:4] - ECA.2)/ECA.2) * 100
      DECA.2[,c(2:3,5:7)] <- round(DECA.2[,c(2:3,5:7)], 3)
      DECA.3 <- ddply(DECA.2, .(scenary), dplyr::summarize,
                      Type_Change = dECAfun(.data$dECA, .data$dA))
      DECA.3$Type <- ddply(DECA.2, .(scenary), dplyr::summarize,
                           Type = dECAfun2(.data$dECA, .data$dA))[[2]]
      names(DECA.3)[2:3] <- c("dA/dECA comparisons", "Type of change")
      DECA.4 <- cbind(DECA.2, DECA.3)
      DECA.4[8] <- NULL

      names(DECA.4)[c(1:3,5)] <- c("Scenary", paste0("Area (",area_unit,")"),
                                   paste0("ECA (",area_unit,")"), "Normalized ECA")
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
    setwd(ttt.2)
    ###plot

    if(!is.null(plot)){
      ECAplot <- map(ECA3, function(x){
        ECA4 <- (x[2] * 100)/ LA #Porcentaje de Area
        names(ECA4) <- "Habitat"
        ECA4$Loss <- 100 - ECA4$Habitat # Perdido
        ECA4$Connected <- x[[5]]
        ECA4$Year <- plot

        #Table 1
        ECA5 <- ECA4
        ECA4 <- melt(ECA4, id="Year")
        names(ECA4)[3] <- "percentage"
        ECA4$variable <- factor(ECA4$variable, levels = c("Loss", "Habitat", "Connected"))

        #Table 2
        ECA5$Connected <- (ECA5$Connected * ECA5$Habitat)/100
        ECA5$Habitat <- ECA5$Habitat - ECA5$Connected
        ECA5 <- melt(ECA5, id="Year")
        names(ECA5)[3] <- "percentage"
        ECA5$variable <- factor(ECA5$variable, levels = c("Loss", "Habitat", "Connected"))
        ECA5$Text <- ECA4$percentage

        ECA5$pos <- ECA5$percentage/2
        ECA5$pos[which(ECA5$variable == "Habitat")] <- (ECA5$percentage[which(ECA5$variable == "Habitat")]/2) + ECA5$percentage[which(ECA5$variable == "Connected")]
        ECA5$pos[which(ECA5$variable == "Loss")] <- (ECA5$percentage[which(ECA5$variable == "Loss")]/2) + ECA5$Text[which(ECA5$variable == "Habitat")]

    #Plot
   fill <- c("#443A82", "#30678D", "#26AE7F")

    if (distance$type  %in% c("centroid", "edge", "hausdorff edge")){
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
            legend.text = element_text(size = 11),
            legend.spacing.x = unit(0.5, 'cm')) +
      labs(x = "Year", y = "Percentage (%)") +
      ggtitle(paste0("Equivalent Connected Area: Dispersal distance = ", dp_text)) +
      scale_fill_manual(values = fill) +
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
    if(is.null(plot)){
      ECA3 <- list(dECA_table = ECA3[[1]][[1]])
    } else {
      ECA3 <- list(dECA_table = ECA3[[1]][[1]], dECA_Plot = ECA3[[1]][[2]])
    }
  }
    }
  return(ECA3)
  }

