#' Test the ECA or ProtConn metrics using multiple dispersal distances
#'
#' @param nodes Object of class sf, sfc, sfg or SpatialPolygons.
#' @param attribute character. Column name with the nodes attribute. If NULL, then the nodes area will be estimated and used as the attribute.
#' @param distance1 list. Distance parameters. For example: type, resistance,or keep. For "type" choose one of the distances: "centroid" (faster), "edge",
#' "least-cost distance" or "commute distance". If the type is equal to "least-cost distance" or "commute distance", then you have to use the "resistance" argument. "keep" is a numeric value used for higher processing.
#'   To See more options consult the help function of distancefile().
#' @param distance2 list. see distance1 argument
#' @param distance3 list. see distance1 argument
#' @param distance4 list. see distance1 argument
#' @param metric character. "IIC" to estimate the ECA using the IIC index, "PC" to estimate the ECA using the PC index or "ProtConn" to estimate the ProtConn indicator using the PC index.
#' @param probability numeric. numeric. Connection probability to the selected distance threshold, e.g., 0.5
#' (default) that is 50 percentage of probability connection. Use in case of selecting the "PC"
#' metric or "ProtConn". If probability = NULL, then it will be the inverse of the mean dispersal distance
#' for the species (1/α; Hanski and Ovaskainen 2000).
#' @param distance_thresholds numeric. Distances thresholds (minimum 3) to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000, 100000); sequence distances (recommended): distance_thresholds = seq(10000,100000, 10000).
#' @param region object of class sf, sfc, sfg or SpatialPolygons. If metric is equal to "ProtConn" then you must provide a region shapefile.
#' @param transboundary numeric. Buffer to select polygones in a second round, their attribute value = 0, see Saura et al. 2017. You can set one transboundary value or one per each threshold distance.
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to m2).
#' @param groups Selected representative threshold distances (distance just before the biggest changes in connectivity metric)
#' @param area_unit character. If attribute is NULL you can set an area unit, "Makurhini::unit_covert()"
#' compatible unit(e.g., "m2", "km2", "ha"). Default equal to hectares "ha".
#' @param write character. Folder path and prefix, for example: "C:/Folder/test".
#' @param intern logical. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @references Correa Ayram, C. A., Mendoza, M. E., Etter, A., & Pérez Salicrup, D. R. (2017). Anthropogenic impact on habitat connectivity: A multidimensional human footprint index evaluated in a highly biodiverse landscape of Mexico. Ecological Indicators, 72, 895–909. https://doi.org/10.1016/j.ecolind.2016.09.007
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(rgeos)
#'
#' data("list_forest_patches", package = "Makurhini")
#' data("study_area", package = "Makurhini")
#'
#' Max_attribute <- unit_convert(gArea(study_area), "m2", "ha")
#'
#' test_metric_distance(nodes = list_forest_patches[[1]],
#'                      distance1 =list(type= "centroid"),
#'                      distance2 =list(type= "edge"),
#'                      attribute = NULL, area_unit = "ha",
#'                      LA = Max_attribute ,
#'                      distance_thresholds = seq(10000,100000, 10000),
#'                      groups = 0)
#'
#'
#' data("Protected_areas", package = "Makurhini")
#' data("regions", package = "Makurhini")
#' region <- regions[1,]
#'
#' test_metric_distance(nodes = Protected_areas,
#'                      distance1 =list(type= "centroid"),
#'                      distance2 =list(type= "edge", keep = 0.05),
#'                      metric = "ProtConn", probability = 0.5,
#'                      area_unit = "ha",
#'                      region = region, transboundary = 50000,
#'                      distance_thresholds = seq(10000,100000, 10000))
#'}
#' @importFrom magrittr %>%
#' @importFrom progressr handlers handler_pbcol progressor
#' @importFrom crayon bgWhite white bgCyan
#' @importFrom ggplot2 ggplot geom_line geom_point aes theme scale_x_continuous scale_y_continuous element_blank labs scale_colour_manual element_line element_text element_rect
#' @importFrom purrr compact map_dfr
#' @importFrom methods as
#' @importFrom utils tail
#' @importFrom sf st_as_sf st_zm st_buffer st_cast st_intersection st_geometry st_difference st_area
test_metric_distance <- function(nodes,
                              attribute = NULL,
                              distance1 = NULL,
                              distance2 = NULL,
                              distance3 = NULL,
                              distance4 = NULL,
                              metric = "IIC", probability = NULL,
                              distance_thresholds,
                              region = NULL,
                              LA = NULL,
                              transboundary = NULL,
                              area_unit = "ha",
                              groups = 3, write = NULL,
                              intern = TRUE){
  if(isFALSE("ggplot2" %in% rownames(installed.packages())) &
     isFALSE("ggpubr" %in% rownames(installed.packages()))){
    stop("To make the plots you need to install the packages ggplot2 and ggpubr")
  }

  if(class(nodes)[1] == "sf") {
    nodes <- as(nodes, 'Spatial')
  }
  options(warn = -1)
  . = NULL
  distances_test <- list(distance1, distance2, distance3, distance4)
  distances_test <- compact(distances_test)

  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if(metric=="ProtConn"){
    if (missing(region)) {
      stop("error missing file of region")
    } else {
      if (is.numeric(region) | is.character(region)) {
        stop("error missing file of region")
      }
    }

    if (is.null(probability) | !is.numeric(probability)) {
      stop("error missing probability")
    }

    region <- st_as_sf(region) %>% st_zm() %>%
      st_buffer(., dist = 0) %>% st_cast("POLYGON")

    nodes <- st_as_sf(nodes) %>% st_zm() %>%
      st_buffer(., dist = 0) %>% st_cast("POLYGON")

    nodes$IDTemp <- 1:nrow(nodes)
    nodes <- nodes[, "IDTemp"]

    nodes.2 <- tryCatch(Protconn_nodes(x = region,
                                       y = nodes,
                                       buff = transboundary,
                                       method = "nodes",
                                       xsimplify = FALSE,
                                       metrunit = area_unit,
                                       protconn_bound = FALSE,
                                       delta = FALSE), error = function(err)err)
    if(is.null(LA)){
      LA <- sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit))
    }

    nodes <- as(nodes.2[[1]], 'Spatial')
    attribute <- "attribute"
    }

  if (metric=="PC"){
    if (is.null(probability) | !is.numeric(probability)) {
      stop("error missing probability")
    }
  }

  if(is.null(groups)){
    groups <- 0
  }

  #Id
  nodes@data$IdTemp <- 1:nrow(nodes)

  loop <- 1:length(distances_test)

  if(isTRUE(intern)){
    handlers(global = T)
    handlers(handler_pbcol(complete = function(s) crayon::bgYellow(crayon::white(s)),
                           incomplete = function(s) crayon::bgWhite(crayon::black(s)),
                           intrusiveness = 2))
    pb <- progressor(along = loop)
  }

  conn_metric <- map_dfr(loop, function(x){
    x.1 <- distances_test[[x]]
    ECA_metric <-  map_dfr(distance_thresholds, function(y) {
      tab1 <- tryCatch(MK_dPCIIC(nodes = nodes, attribute = attribute,
                                 restoration = NULL,
                                 distance = x.1, area_unit = area_unit,
                                 metric = if(metric=="ProtConn"){"PC"}else{metric},
                                 probability = probability,
                                 distance_thresholds = y,
                                 overall = TRUE, onlyoverall = TRUE,
                                 LA = LA, rasterparallel = FALSE, write = NULL), error = function(err)err)

      if (inherits(tab1, "error")){
        stop(tab1)
      }

      tab1 <- tab1[2,2]

      if(metric=="ProtConn"){
        tab1 <- 100 * (tab1 / LA)
      }
      return(tab1)
    })

    if(length(distances_test)>1 & isTRUE(intern)){
      pb()
    }
    conn_metric2 <- cbind(ECA_metric, distance_thresholds)
    conn_metric2 <- as.data.frame(conn_metric2)
    names(conn_metric2) <- c(if(metric=="ProtConn"){"ProtConn"}else{"ECA"}, "Distance")
    conn_metric2$Group <- x.1$type
    return(conn_metric2)
  })


  ###Differences
  peaks_groups <- list()
  unique_groups <- unique(conn_metric$Group)


  ###Groups return
  if(groups>0){
    for (i in 1:length(unique_groups)){
      table1 <- conn_metric[conn_metric$Group == unique_groups[i],]
      dif1 <- table1[,1]

      change <- list()

      for(j in 2:length(dif1)){
        change[[j]] <- (abs(dif1[j]-dif1[j-1])*100)/dif1[j]
      }
      change <- do.call(rbind, change) %>% as.data.frame()
      change$metric <- table1[2:nrow(table1),1]
      change$distance <- table1$Distance[2:nrow(table1)]
      change <- change[order(change$V1),]
      change <- change[(nrow(change)-(groups-1)):nrow(change),]
      change$Group <- 1:groups
      change <- change[order(change$distance),]

      peak <- data.frame(Group = change$Group)
      peak$Dist_from <- c(min(table1$Distance), format(change$distance[1:(groups-1)], scientific = F))
      peak$Dist_to <- c(format(change$distance[1:(groups)], scientific = F))
      peak$Dist_to[groups] <- format(max(table1$Distance),  scientific = F)
      peak$metric <- change$metric

      names(peak)[4] <- if(metric=="ProtConn"){"ProtConn"}else{"ECA"}
      if(length(unique(peak[,4])) == 1){
        peak <- peak[which(peak$Group ==1),]
      }
      peaks_groups[[i]] <- peak
    }

    names(peaks_groups) <- unique_groups

    ###Groups plot
    peak_plot <- list()
    for(i in 1:length(unique_groups)){
      peaks_p <- cbind(peaks_groups[[i]], unique_groups[[i]])
      names(peaks_p)[c(3, 5)] <- c("To", "plot_group")
      peak_plot[[i]] <- peaks_p
    }
    peak_plot <- do.call(rbind, peak_plot)
  }
  ###Plot
  if(is.null(attribute)){
    ytitle = if(metric=="ProtConn"){"ProtConn (%)"}else{"ECA (ha)"}
  } else {
    ytitle = if(metric=="ProtConn"){"ProtConn (%)"}else{paste0("ECA", " (", area_unit, ")")}
  }


  ccolour <- c("#E16A86", "#909800", "#00AD9A", "#9183E6")
  ccolour <- ccolour[1:length(unique_groups)]

  p1 <- ggplot() +
    geom_line(aes(y = conn_metric[,1], x = conn_metric$Distance, colour = conn_metric$Group),
              size = 1.5, data = conn_metric,
              stat="identity") +
    theme(legend.position="bottom", legend.direction="horizontal",
          legend.title = element_blank()) +
    scale_x_continuous(breaks = c(seq(min(conn_metric$Distance), max(conn_metric$Distance),
                                      max(conn_metric$Distance)/10)[2:10] -min(conn_metric$Distance),
                                  max(conn_metric$Distance)),
                       labels = function(x) format(x, scientific = FALSE)) +
    scale_y_continuous(breaks = c(seq(min(conn_metric[,1]), max(conn_metric[,1]),
                                      (max(conn_metric[,1])-min(conn_metric[,1]))/4)),
                       labels = function(x) round(x, 2))+
    labs(x = "Dispersal distance", y = ytitle) +
    scale_colour_manual(values = ccolour) +
    theme(axis.line = element_line(size = 1, colour = "black"),
          panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) +
    theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          text = element_text(family = "Tahoma", size = 12),
          axis.text.x= element_text(colour = "black", size = 10, angle = 90),
          axis.text.y= element_text(colour = "black", size = 10),
          legend.key= element_rect(fill = "white", colour = "white"),
          legend.text = element_text(size = 11))

  if(groups>0){
    peak_plot$To <- as.numeric(peak_plot$To)
    p1 <- p1 + geom_point(data = peak_plot,
                          mapping = aes(x = conn_metric[which(conn_metric[,1] %in% peak_plot[,4]),2],
                                        y = peak_plot[,4]), size = 4)+
      geom_text(aes(x = conn_metric[which(conn_metric[,1] %in% peak_plot[,4]),2],
                    y = peak_plot[,4],
                    label = conn_metric[which(conn_metric[,1] %in% peak_plot[,4]),2]),
                data = peak_plot, hjust = -0.3, vjust = -0.5)+
      theme(legend.text = element_text(size = 11))
  }

  result <- list(conn_metric,if(groups>0){peaks_groups}, p1)
  result <- purrr::compact(result)

  if(groups>0){
    names(result) <- c("Table_metric_test",  "Table_groups_test", "plot_metric_test")
  }else {
    names(result) <- c("Table_metric_test", "plot_metric_test")
  }

  if(!is.null(write)){
    ggsave(paste0(write, '_metricTest.tif'), plot = p1, device = "tiff", width = 12,
           height = 8, compression = "lzw")
    names(peak_plot)[3] <- "To(dispersal distance)"

    write.csv(conn_metric, paste0(write, '_MetricValTest.csv'), row.names = FALSE)

    if(groups>0){
      write.csv(peak_plot, paste0(write, '_GroupsTest.csv'), row.names = FALSE)
    }
  }
  return(result)
}





