#' Fragmentation Statistics
#'
#' Calculate patches/nodes and landscape statistics
#' @param nodes Object of class \code{sf, sfc, sfg, SpatialPolygons}. Individual nodes, the object must be in a projected coordinate system.
#' @param edge_distance \code{numeric}. Distance to edge in meters. Default equal 500 m (Haddad et al. 2015)
#' @param min_node_area \code{numeric}. Minimum node area used to calculate the number of nodes with an area smaller than the one provided. Default equal 100 km\out{<sup>2</sup>} (Haddad et al. 2015)
#' @param landscape_area \code{numeric}. Total area of the study landscape in km\out{<sup>2</sup>} (optional). If NULL the total nod area will be used.
#' @param area_unit \code{character}. You can set an area unit (e.g., "km2", "cm2", "m2", "ha"; see \link[Makurhini]{unit_convert}). Default equal to square kilometers "km2".
#' @param perimeter_unit \code{character}. You can set a perimeter unit (e.g., "km", "cm", "m", "ha"; see \link[Makurhini]{unit_convert}). Default equal to kilometers "km".
#' @param plot \code{logical}. Basic histograms and core area - edge map.
#' @param write \code{character}. Write the table (landscape statistics), sf object (patch/node statistics) and plots. It's necessary to specify the path and prefix, for example,
#' to save in the path "C:/Folder" with the prefix "Fragmentation": \code{"C:/Folder/Fragmentation"}.
#' @return
#' Patch/node and landscape statistics:\cr
#' 1) Patches Area in square kilometers.\cr
#' 2) Number patches.\cr
#' 3) Mean size of patches.\cr
#' 4) Number of patches smaller than the parameter \code{min_node_area} (default = 100 km\out{<sup>2</sup>}).\cr
#' 5) Percentage of patches smaller than the parameter \code{min_node_area} (default = 100 km\out{<sup>2</sup>}).\cr
#' 6) Total patch edge. Total perimeter of the patches (unit = \code{perimeter_unit}).\cr
#' 7) Edge density. Total perimeter per unit of area (unit = \code{area_unit}), default = km\out{<sup>2</sup>}. A value of 0 is present when there is no edge in the landscape.\cr
#' 8) Patch density.\cr
#' 9) Total core area (units = \code{area_unit}) considering the distance set in the parameter \code{edge_distance} (delfault = 500 m).\cr
#' 10) Cority index. It is a measure of fragmentation with respect to a distance from the core area (parameter \code{edge_distance}; delfault = 500 m), where a value of 1 indicates a landscape without fragmentation. Average for landscape level.\cr
#' 11) Shape Index. A simple shape metric that takes values from 1 (perfectly compact) to infinity is derived by dividing the perimeter by the perimeter of a circle of the same area. Average for landscape level. \cr
#' 12) Fractal dimension. The index reflects the complexity of the shape of the fragment. A fractal dimension greater than 1 indicates an increase in the complexity of the shape. When the value is close to 1 the shape is simple, such as squares.Average for landscape level. \cr
#' 13) Effective Mesh Size. Effective Mesh Size (MESH) is a measure of the degree of fragmentation in the landscape ranging from 0 to the total landscape area. MESH is maximum when the landscape unit consists of a single habitat fragment or the habitat is continuous beyond the landscape unit analyzed (Moser, 2007).\cr
#' 14) Core percent (patch level). Percentage of core area in the patch (units = \code{area_unit}) considering the distance set in the parameter \code{edge_distance} (delfault = 500 m).\cr
#' 15) Edge percent (patch level). Percentage of edge in the patch (units = \code{area_unit}) considering the distance set in the parameter \code{edge_distance} (delfault = 500 m).\cr
#' 16) PARA (patch level). Ratio of the patch perimeter to area.
#' *NOTE.* In the results we use the term patches instead of nodes due to the common use of this term in fragmentation statistics in science.
#' @references
#' Haddad et al. (2015). Science Advances 1(2):e1500052. https://www.science.org/doi/10.1126/sciadv.1500052.\cr
#' McGarigal, K., S. A. Cushman, M. C. Neel, and E. Ene. 2002. FRAGSTATS: Spatial Pattern Analysis Program for Categorical Maps. Computer software program produced by the authors at the University of Massachusetts, Amherst. Available at the following web site:
#'  \url{www.umass.edu/landeco/research/fragstats/fragstats.html}.\cr
#' Moser, B., Jaeger, J.A.G., Tappeiner, U. et al. Modification of the effective mesh size for measuring landscape fragmentation to solve the boundary problem. Landscape Ecol 22, 447–459 (2007).  \url{https://doi.org/10.1007/s10980-006-9023-0}
#' @examples
#' data("habitat_nodes", package = "Makurhini")
#' nrow(habitat_nodes) # Number of nodes
#' fragmentation <- MK_Fragmentation(nodes = habitat_nodes, edge_distance = 1000, plot = TRUE)
#' #Table
#' fragmentation$`Summary landscape metrics (Viewer Panel)`
#' #Shapefile
#' fragmentation$`Patch statistics shapefile`
#' @export
#' @importFrom magrittr %>%
#' @importFrom sf st_as_sf st_zm st_cast st_buffer st_area st_length st_boundary st_geometry st_geometry<- write_sf st_is_empty
#' @importFrom raster plot
#' @importFrom ggplot2 ggplot aes geom_histogram theme element_blank element_text labs ggsave geom_sf theme_bw scale_fill_manual
#' @importFrom formattable formattable formatter style
#' @importFrom ggpubr ggarrange
#' @importFrom utils installed.packages
MK_Fragmentation <- function(nodes = NULL, edge_distance = 500, min_node_area = 100,
                             landscape_area = NULL, area_unit = "ha", perimeter_unit = "km",
                             plot = FALSE, write = NULL){
  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
    } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(nodes)[1] == "SpatialPolygonsDataFrame") {
    nodes <- st_as_sf(nodes) %>% st_cast("POLYGON")
  }

  colors <- c("Core" = "#1a9641", "Edge" = "Red")
  vcol <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")

  nodes$IdTemp <- 1:nrow(nodes)

  ###Patch metrics
  CoreA <- st_buffer(nodes, dist = -(edge_distance))
  data <- data.frame(IdTemp = nodes$IdTemp,
                     Area = round(unit_convert(st_area(nodes, by_element = TRUE), "m2", area_unit), 4),
                     CA = round(unit_convert(st_area(CoreA), "m2", area_unit), 4))
  data$CAPercent <- round((data$CA * 100) / data$Area, 4)
  data$Perimeter <- round(unit_convert(st_length(st_boundary(nodes)), "m", perimeter_unit), 3) %>% as.numeric()
  data$EdgePercent <- round((100 - data$CAPercent), 4)
  data$PARA <- round(data$Area / data$Perimeter, 4)
  data$ShapeIndex <- round(data$Perimeter / (2 * sqrt(data$Area*pi)), 4)
  data$FRAC <- round((2 * (log(data$Perimeter))) / log(data$Area), 4)
  patches <- base::merge(nodes, data, by = "IdTemp", all = T)
  patches$IdTemp <- NULL

  #Landscape metrics
  if (is.null(landscape_area)){
    landscape_area <- round(sum(data$Area, na.rm = TRUE), 4)
  }

  LM <- data.frame(a = round(sum(data$Area, na.rm = TRUE), 4),
                 b = nrow(data),
                 c = round(mean.default(data$Area, na.rm = TRUE), 4),
                 d = length(which(data$Area < min_node_area)),
                 e = round(sum(data[which(data$Area < min_node_area),]$Area, na.rm = TRUE) * 100 / sum(data$Area, na.rm = TRUE), 4),
                 f = round(sum(data$Perimeter, na.rm = TRUE), 4),
                 g = round(sum(data$Perimeter, na.rm = TRUE) / sum(data$Area, na.rm = TRUE), 4),
                 l = round(((nrow(data)/unit_convert(landscape_area, area_unit, "m2")) * unit_convert(1, area_unit, "m2"))*100, 4),
                 h = round(sum(data$CA, na.rm = TRUE), 4),
                 i = round((nrow(data) - length(which(data$CA == 0))) / nrow(data), 4),
                 j = round(mean(data$ShapeIndex, na.rm = TRUE), 4),
                 k = round(mean(data$FRAC, na.rm = TRUE), 4),
                 m = round((1/landscape_area) * sum(data$Area^2), 4))

  LM_names <- c(paste0("Patch area (", area_unit, ")"),
                "Number of patches",
                "Size (mean)",
                "Patches < minimum patch area",
                "Patches < minimum patch area (%)",
                "Total edge", "Edge density", "Patch density",
                paste0("Total Core Area (", area_unit, ")"), "Cority",
                "Shape Index (mean)",
                "FRAC (mean)", paste0("MESH (", area_unit, ")"))
  LM <- t(LM) %>% as.data.frame(); LM$Metric <- LM_names
  LM$Value <- LM[[1]]; rownames(LM) <- NULL; LM[1] <- NULL

  #Plot

  if(isTRUE("ggplot2" %in% rownames(installed.packages())) &
     isTRUE("ggpubr" %in% rownames(installed.packages()))){
    data2 <- data
    data2$Area <- ifelse(data$Area <= 0, 0, log10(data$Area))
    data2$CA <- ifelse(data$CA <= 0, 0, log10(data$CA))
    data2$Perimeter <- ifelse(data$Perimeter <= 0, 0, log10(data$Perimeter))
    if(!is.null(write) & isTRUE(plot)) {
      par(mfrow = c(1,1))
      p0 <- ggplot(data = patches) +
        geom_sf(colour = "Red", aes(fill = "Edge"))+
        geom_sf(data = CoreA[which(!st_is_empty(CoreA)),], colour = "#1a9641", aes(fill = "Core"))+
        scale_fill_manual(name = "Legend", values = colors)+
        theme_bw(base_size = 22)

      ggsave(paste0(write, '_CoreEdge.tif'), plot = p0, device = "tiff", width = 15,
             height = 11, compression = "lzw")

      p1 <- ggplot(data2, aes(x = data2$Area)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = paste0("log10 (", area_unit, ")"), y = "Frequency", title = "Size") +
        theme(plot.title = element_text(size=20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))+
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p2 <- ggplot(data2, aes(x = data2$Perimeter)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T)+
        labs(x = paste0("log10 (", perimeter_unit, ")"), y ="Frequency", title = "Perimeter") +
        theme(plot.title = element_text(size = 20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))+
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p3 <- ggplot(data2, aes(x = data2$ShapeIndex)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "Shape Index", y = "Frequency", title = "Shape Index") +
        theme(plot.title = element_text(size = 20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))+
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p4 <- ggplot(data2, aes(x = data2$CA)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T)+
        labs(x = paste0("log10 (", area_unit, ")"), y = "Frequency", title = "Core Area") +
        theme(plot.title = element_text(size = 20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20)) +
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p5 <- suppressWarnings(ggarrange(p1, p2, p3, p4))
      ggsave(paste0(write, '_fragStats.tif'),  plot = p5, device = "tiff", width = 15,
             height = 11, compression = "lzw")
    }
    if (isTRUE(plot)) {
      p0 <- ggplot(data = patches) +
        geom_sf(colour = "Red", aes(fill = "Edge"))+
        geom_sf(data = CoreA[which(!st_is_empty(CoreA)),], colour = "#1a9641", aes(fill = "Core"))+
        scale_fill_manual(name = "Legend", values = colors)+
        theme_bw()
      p1 <- ggplot(data2, aes(x = data2$Area)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = paste0("log10 (", area_unit, ")"), y = "Frequency", title = "Size") +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      p2 <- ggplot(data2, aes(x = data2$Perimeter)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = paste0("log10 (", perimeter_unit, ")"), y ="Frequency", title = "Perimeter") +
        theme(plot.title = element_text(size=14, face="bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      p3 <- ggplot(data2, aes(x = data2$ShapeIndex)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "Shape Index", y ="Frequency", title = "Shape Index") +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      p4 <- ggplot(data2, aes(x = data2$CA)) +
        geom_histogram(color = "black", fill = vcol,
                       bins = 10, position = "dodge", na.rm = T) +
        labs(x = paste0("log10 (", area_unit, ")"), y ="Frequency", title = "Core Area")+
        theme(plot.title = element_text(size=14, face ="bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))

      p5 <- suppressWarnings(ggarrange(p1, p2, p3, p4))

    }
    rm(data2)
  } else {
    message("To make the plots you need to install the packages ggplot2 and ggpubr")
    plot = FALSE
  }

  #Write outputs
  if (!is.null(write)) {
    write.csv(LM, paste0(write, '_LandscapeMetrics.csv'))
    write.csv(data, paste0(write, '_PatchMetrics.csv'))
    n <- which(names(patches) %in% c("Area", "CA", "CAPercent", "Perimeter", "EdgePercent", "PARA",
                          "ShapeIndex", "FRAC"))
    patches_2 <- patches
    names(patches_2)[n] <- c("Area", "CA", "CA%", "P", "Edge%", "PARA", "ShapeI", "FRAC")
    write_sf(patches_2, paste0(write, '_PatchMetrics.shp'), delete_layer = TRUE)
  }

  #Return
  LM <- formattable(LM,
                    align = c("l","c"),
                    list(`Indicator Name` = formatter("span",
                                                      style = ~ style(color = "grey", font.weight = "bold"))))

  if (isTRUE(plot)){
    base::print(p0); base::print(p5)
    return(list("Summary landscape metrics (Viewer Panel)" = LM,
                "Patch statistics shapefile" = patches))
  } else {
    return(list("Summary landscape metrics (Viewer Panel)" = LM,
                "Patch statistics shapefile" = patches))
    }
  }

