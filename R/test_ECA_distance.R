#'Test the ECA metric to multiple dispersal distances
#'
#' @param nodes Object of class sf, sfc, sfg or SpatialPolygons.
#' @param attribute character. Column name with the nodes attribute. If NULL, then the nodes area will be estimated and used as the attribute.
#' @param distance1 list. Distance parameters. For example: type, resistance,or keep. For "type" choose one of the distances: "centroid" (faster), "edge",
#' "least-cost distance" or "commute distance". If the type is equal to "least-cost distance" or "commute distance", then you have to use the "resistance" argument. "keep" is a numeric value used for higher processing.
#'   To See more options consult the help function of distancefile().
#' @param distance2 list. see distance1 argument
#' @param distance3 list. see distance1 argument
#' @param distance4 list. see distance1 argument
#' @param metric character. "IIC" considering topologycal distances or "PC" considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "BPC" metric.
#' @param distance_thresholds numeric. Distances thresholds (minimum 3) to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000, 100000); sequence distances (recommended): distance_thresholds = seq(10000,100000, 10000).
#' @param groups Selected representative threshold distances (distance just before the biggest changes in connectivity metric)
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to m2).
#' @param write character. Write output shapefile, example, "C:/ejemplo.shp".
#' @references Correa Ayram, C. A., Mendoza, M. E., Etter, A., & Pérez Salicrup, D. R. (2017). Anthropogenic impact on habitat connectivity: A multidimensional human footprint index evaluated in a highly biodiverse landscape of Mexico. Ecological Indicators, 72, 895–909. https://doi.org/10.1016/j.ecolind.2016.09.007
#' @export
#' @examples
#' \dontrun{
#' ruta <- system.file("extdata", "ECA_example.RData", package = "Makurhini")
#' load(ruta)
#' test_ECA_distance(nodes = forest_patches[[1]], distance1 =list(type= "centroid"), LA = 279165,
#'                   distance_thresholds =  seq(10000,100000, 10000))
#'                   }
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom  purrr compact
#' @importFrom methods as
#' @importFrom utils tail
test_ECA_distance <- function(nodes,
                              attribute = NULL,
                              distance1 = NULL,
                              distance2 = NULL,
                              distance3 = NULL,
                              distance4 = NULL,
                              metric = "IIC", probability = NULL,
                              distance_thresholds, LA = NULL,
                              groups = 3, write = NULL){
  if(class(nodes)[1] == "sf") {
    nodes <- as(nodes, 'Spatial')
    }
  options(warn = -1)
  ttt.2 <- getwd()
  #
  temp.1 <- paste0(tempdir(), "/TempInputs", sample(1:1000, 1, replace = TRUE))
  dir.create(temp.1, recursive = T)

  #Id
  nodes@data$IdTemp <- 1:nrow(nodes)

  #Nodes
  if(is.null(attribute)){
    nodesfile(nodes, id = "IdTemp", attribute = NULL, area_unit = "ha", write = paste0(temp.1,"/nodes.txt"))
  } else {
    nodesfile(nodes, id = "IdTemp", attribute = attribute, write = paste0(temp.1,"/nodes.txt"))
  }

  setwd(temp.1)

  distances_test <- list(distance1, distance2, distance3, distance4)
  distances_test <- compact(distances_test)

  conn_metric <- lapply(distances_test, function(x){

    distancefile(nodes,  id = "IdTemp", type = x$type, keep = x$keep,
                 resistance = x$resistance, CostFun = x$CostFun, ngh = x$ngh,
                 threshold = x$threshold,
                 distance_unit = x$distance_unit, x$geometry_out,
                 write = paste0(temp.1,"/Dist.txt"))

    if (is.null(x$threshold)) {
      pairs = "all"
    } else {
      pairs = "notall"
    }

    conn_metric2 <- lapply(as.list(distance_thresholds), function(y){
      tab1 <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                         typeconnection = "dist", typepairs = pairs,
                         index = metric, thdist = y,
                         multdist = NULL, conprob = probability,
                         onlyoverall = TRUE, LA = LA,
                         nrestauration = FALSE,
                         prefix = NULL, write = NULL)
      tab1 <- tab1[[which(lapply(tab1, function(x) paste0(nrow(x), ncol(x))) == "32" | lapply(tab1, function(x) paste0(nrow(x), ncol(x))) == "22")]]
      tab1 <- tab1[[2,2]]
      return(tab1)})

    conn_metric2 <- do.call(rbind, conn_metric2)
    conn_metric2 <- cbind(conn_metric2, distance_thresholds)
    conn_metric2 <- as.data.frame(conn_metric2)
    names(conn_metric2) <- c("ECA", "Distance")
    conn_metric2$Group <- x$type
    return(conn_metric2)
    })

  setwd(ttt.2)
  #Join ECA tables
  conn_metric <-  do.call(rbind, conn_metric)

  #Selection
  ###Differences
  peaks_groups <- list()
  unique_groups <- unique(conn_metric$Group)

  ###Groups return
  for (i in 1:length(unique_groups)){
      table1 <- conn_metric[conn_metric$Group == unique_groups[i],]
      dif1 <- table1$ECA[2:length(table1$ECA)]
      dif2 <- ((dif1 * 100)/table1$ECA) - 100
      dif2 <-  round(dif2, 6)
      peak <- table1[which(dif2 %in% tail(sort(dif2),(groups - 1))) +1,]
      From <- c(min(table1$Distance), peak[[2]] + 1)
      To <- c(peak[[2]], max(table1$Distance))
      group_i <- 1:groups
      group_i <- as.data.frame(cbind(group_i, From, To))
      names(group_i)[c(1, 3)] <- c("Group", "To(dispersal distance)")
      group_i$ECA <- c(peak$ECA, table1[which(table1$Distance == max(table1$Distance)), 1])
      peaks_groups[[i]] <- group_i
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

    ###Plot
    if(is.null(attribute)){
      ytitle = "ECA (ha)"
      } else {
        ytitle = "ECA"
      }


    ccolour <- c("#E16A86", "#909800", "#00AD9A", "#9183E6")
    ccolour <- ccolour[1:length(unique_groups)]

    p1 <- ggplot() +
        geom_line(aes(y = conn_metric$ECA, x = conn_metric$Distance, colour = conn_metric$Group),
                  size = 1.5, data = conn_metric,
                  stat="identity") +
        theme(legend.position="bottom", legend.direction="horizontal",
              legend.title = element_blank()) +
        scale_x_continuous(breaks = c(seq(min(conn_metric$Distance), max(conn_metric$Distance),
                                          max(conn_metric$Distance)/10)[2:10] -min(conn_metric$Distance),
                                      max(conn_metric$Distance))) +
        scale_y_continuous(breaks = c(seq(min(conn_metric$ECA), max(conn_metric$ECA),
                                          (max(conn_metric$ECA)-min(conn_metric$ECA))/4)))+
        labs(x = "Dispersal distance", y = ytitle) +
        scale_colour_manual(values = ccolour) +
        theme(axis.line = element_line(size = 1, colour = "black"),
              panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank()) +
        theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
              text = element_text(family = "Tahoma", size = 12),
              axis.text.x= element_text(colour = "black", size = 10, angle = 90),
              axis.text.y= element_text(colour = "black", size = 10),
              legend.key= element_rect(fill = "white", colour = "white"))+
        geom_point(data = peak_plot,
                   mapping = aes(x = peak_plot$To, y = peak_plot$ECA), size = 4)+
        geom_text(aes(x = peak_plot$To, y = peak_plot$ECA, label = peak_plot$To),
                  data = peak_plot, hjust = -0.3, vjust = -0.5)+
        theme(legend.text = element_text(size = 11))

    result <- list(peaks_groups, p1)
    names(result) <- c("Table_ECA_test", "plot_ECA_test")

    if(!is.null(write)){
      ggsave(paste0(write, '_ECATest.tif'), plot = p1, device = "tiff", width = 12,
             height = 8, compression = "lzw")
      names(peak_plot)[3] <- "To(dispersal distance)"
      write.table(peak_plot, paste0(write, '_ECATest.txt'),
                  sep = "\t", row.names = FALSE)
    }
   return(result)
}

#'Plot probability  of dispersal
#'
#' Negative exponential dispersal kernel to calculate the probability of dispersal between two nodes. Used in probabilistic indexes (e.g., PC, BCPC, ProtConn etc.)
#' @param probability numeric. Probability of dispersal at max_distance. If NULL
#' @param max_distance numeric. Up to five maximum dispersal distance (meters)
#' @param eval_distance numeric. Calculate the probability of dispersal at a specific distance
#' @param min.prob numeric. Value between 0-1, maximum x axe value
#' @examples
#' \dontrun{
#' probability_distance(probability= 0.5, max_distance = c(1000, 10000, 30000, 100000))
#' probability_distance(probability= 0.5, max_distance = 30000, eval_distance = 10000)
#' }
#' @export
#' @importFrom purrr map
#' @importFrom graphics par plot axis box lines legend

probability_distance <- function(probability, max_distance, eval_distance = NULL, min.prob = 0.2){
  x <-1
  repeat{
    x = x+1
    m <- exp(x * log(probability)/max(max_distance))
    if (m < min.prob & x > max(max_distance)){
      break
      return(x)
    }
  }

  m = x
  Rcolors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

  if(!is.null(eval_distance)){
    if(m < eval_distance){
      m <- eval_distance + eval_distance/100
      }
    if(length(max_distance) == 1){
      data1 <- exp(1:m * log(probability)/max_distance)

      result_1 <- exp(eval_distance * log(probability)/max_distance)

      t <- c(rep(result_1, eval_distance - 1), result_1, 0)

      par(xaxs="i")
      plot(data1, type = "l", xlab= "Distance (dij)", ylab = "Probability of dispersal (pij)",
           lwd = 2, xlim=c(0,m), axes = F)
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()
      lines(t, type = "l", col = Rcolors[1], lty = 2, lwd = 2)
      legend("topright", c(paste0("d ", max_distance), "eval_distance"), lwd=c(2,2), lty = c(1,2), col=c("black",Rcolors[1]), y.intersp=1.5)
      return(result_1)
    } else {
      data_list <- map(as.list(max_distance), function(x){
        data1 <- exp(1:m * log(probability)/x)
        result_1 <- exp(eval_distance * log(probability)/x)
        t <- c(rep(result_1, eval_distance - 1), result_1, 0)
        return(list(data1, result_1, t))})
      par(xaxs="i")
      plot(data_list[[1]][[1]], type = "l", xlab= "Distance (dij)", ylab = "Probability of dispersal (pij)",
           lwd = 2, xlim=c(0,m), axes=F)
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()
      lines(data_list[[1]][[3]], type = "l", col = Rcolors[1], lty = 2, lwd = 2)

      for(i in 2:length(max_distance)){
        lines(data_list[[i]][[1]], type = "l", col = Rcolors[i], lwd = 2)
        lines(data_list[[i]][[3]], type = "l", col = Rcolors[1], lty = 2, lwd = 2)
      }
      legend("topright", c(paste0("d ", max_distance), "eval_distance"), lwd=2, lty = c(rep(1,length(max_distance)), 2),
             col= c("black", Rcolors[2:length(max_distance)], Rcolors[1]), y.intersp=1.5)

      result_1 <- map(data_list, function(x){x[[2]]})
      names(result_1) <- paste0("d_", max_distance)
      return(result_1)
    }
    } else {
    if(length(max_distance) == 1){
      data <- exp(1:m * log(probability)/max_distance)
      par(xaxs="i")
      plot(data, type = "l", xlab=  "Distance (dij)", ylab = "Probability of dispersal (pij)",  axes=F)
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()
      lines(data_list[[1]][[3]], type = "l", col = Rcolors[1], lty = 2, lwd = 2)
      legend("topright", paste0("d ", max_distance), lwd=c(2,2), lty = c(1,2), col=c("black",Rcolors[1]), y.intersp=1.5)
      } else {
        data_list <- map(as.list(max_distance), function(x){
          data1 <- exp(1:m * log(probability)/x)
          return(data1)})

        par(xaxs="i")
        plot(data_list[[1]], type = "l", xlab= "Distance (dij)", ylab = "Probability of dispersal (pij)",
             lwd = 2, xlim=c(0,m), axes=F)
        axis(side = 1, at= round(seq(0,m, m/10)))
        axis(side = 2, at= seq(0, 1, 1/10))
        box()

        for(i in 2:length(max_distance)){
          lines(data_list[[i]], type = "l", col = Rcolors[i-1], lwd = 2)
        }
        legend("topright", paste0("d ", max_distance), lwd=2,
               col= c("black", Rcolors[1:length(max_distance)]), y.intersp=1.5)
      }
    }

}





