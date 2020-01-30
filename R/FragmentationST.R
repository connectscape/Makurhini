#' Fragmentation Statistics
#'
#' Calculate patch and landscape statistics
#' @param patches Object of class sf, sfc, sfg or SpatialPolygons with individual patches. The shapefile must be in a projected coordinate system.
#' @param edge_distance Numeric. Distance to edge. Default equal 500 m (Haddad et al. 2015)
#' @param min_patch_area Numeric. Minimum patch area in km². Default equal 100 km²(Haddad et al. 2015)
#' @param landscape_area Numeric. Total landscape area in km² (optional).If NULL the total patch area will be used.
#' @param plot Logical. Basic histograms and core area - edge map.
#' @param write Character. Write the following outputs: Fragmentation.csv, Fragmentation.shp and plots. It's necessary to specify the path and prefix, for example: "C:/Folder/Fragmentation".
#' @return
#' Patch and landscape statistics:\cr
#' 1) Patches Area in square kilometers.\cr
#' 2) Number patches.\cr
#' 3) Mean size of patches.\cr
#' 4) Number of patches smaller than 100 km²or another minimum patch area.\cr
#' 5) Area in square kilometers of patches smaller than 100 km² or another minimum patch area.\cr
#' 6) Total edge.\cr
#' 7) Edge density.\cr
#' 8) Total core area (km²; considering a distance to edge of 500 m).\cr
#' 9) Total core area.\cr
#' 10) Core percent.\cr
#' 11) Edge percent.\cr
#' 12) Cority index. It is a measure of fragmentation with respect to a distance from the core area of 500 m (Haddad et al., 2015), where a value of 1 indicates a landscape without fragmentation.\cr
#' 13) Shape Index. A simple shape metric that takes values from 1 (perfectly compact) to infinity is derived by dividing the perimeter by the perimeter of a circle of the same area.\cr
#' 14) Fractal dimension.\cr
#' 15) Effective Mesh Size.\cr
#' @references McGarigal, K., S. A. Cushman, M. C. Neel, and E. Ene. 2002. FRAGSTATS: Spatial Pattern Analysis Program for Categorical Maps. Computer software program produced by the authors at the University of Massachusetts, Amherst. Available at the following web site:
#'  \url{www.umass.edu/landeco/research/fragstats/fragstats.html}.\cr
#'  \cr
#'  Haddad et al. (2015). Science Advances 1(2):e1500052. \url{DOI: 10.1126/sciadv.1500052}.
#' @examples ruta <- system.file("extdata", "Habitat_Patches.shp", package = "Makurhini")
#'  cores <- sf::read_sf(ruta)
#'  fragmentation <- FragmentationST(patches = cores, edge_distance = 500, plot = TRUE)
#'  names(fragmentation)
#'  #Table
#'  fragmentation$`Summary landscape metrics (Viewer Panel)`
#'  #Shapefile
#'  fragmentation$`Patch statistics shapefile`
#'  #Plot
#'  fragmentation$`Frequency landscape metric plot`
#' @export
#' @importFrom magrittr %>%
FragmentationST <- function(patches, edge_distance = 500, min_patch_area = 100, landscape_area = NULL,
                            plot = FALSE, write = NULL){
  if (missing(patches)) {
    stop("error missing shapefile file of patches")
  } else {
    if (is.numeric(patches) | is.character(patches)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(patches)[1] == "SpatialPolygonsDataFrame") {
    patches <- sf::st_as_sf(patches) %>% sf::st_cast("POLYGON")
  }

  patches$IdTemp <- 1:nrow(patches)

  data <- patches[,which(colnames(patches) == "IdTemp")]
  sf::st_geometry(data) <- NULL

  ###Patch metrics
  CoreA <- sf::st_buffer(patches, dist = -(edge_distance))
  data <- cbind(data, data.frame(Area = cbind(round((sf::st_area(patches, byid = T) * 1e-6), 3)),
              CA = cbind(round(as.numeric(sf::st_area(CoreA) * 1e-6), 3))))
  data$CAPercent <- round((data$CA * 100) / data$Area, 3)
  data$Perimeter <- cbind(round((sf::st_length(sf::st_boundary(patches)) * 0.001), 3))
  data$EdgePercent <- round((100 - data$CAPercent), 3)
  data$PARA <- round(data$Area / data$Perimeter, 3)
  data$ShapeIndex <- round((data$Perimeter / (2 * pi * sqrt(data$Area/pi))), 3)
  data$FRAC <- round((2 * (log(0.25 * data$Perimeter))) / log(data$Area), 3)
  patches <- base::merge(patches, data, by = "IdTemp", all = T)
  patches$IdTemp <- NULL

  ###Landscape metrics
  LM <- data.frame(a = round(sum(data$Area), 3),
                 b = nrow(data),
                 c = round(mean.default(data$Area), 3),
                 d = length(which(data$Area < min_patch_area)),
                 e = round(sum(data[which(data$Area < min_patch_area),]$Area) * 100 / sum(data$Area), 3),
                 f = round(sum(data$Perimeter), 3),
                 g = round(sum(data$Perimeter) / sum(data$Area), 3),
                 h = round(sum(data$CA), 3),
                 i = round((nrow(data) - sum(data$CA == 0, na.rm = TRUE)) /
                                  (sum(data$CA > 0, na.rm = TRUE) +
                                     sum(data$CA == 0, na.rm = TRUE)), 3),
                 j = round(mean(data$ShapeIndex), 3),
                 k = round(mean(data$FRAC), 3))

  LM_names <- c("Patch area (km²)", "Number of patches", "Size (mean)",
                   "Patches < minimum patch area", "Patches < minimum patch area (%)",
                   "Total edge", "Edge density",
                   "Total Core Area (km²)", "Cority",
                   "Shape Index (mean)",
                   "FRAC (mean)", "MESH (km²)")

  if (is.null(landscape_area)){
    Mesh <- data.frame(Mesh = round((1 / LM[1]) * sum(data$Area^2), 3))
    } else {
      Mesh <- data.frame(Mesh = round((1/landscape_area) * sum(data$Area^2), 3))
      }

  LM <- t(cbind(LM, Mesh)) %>% as.data.frame()
  LM$Metric <- LM_names
  LM$Value <- LM[[1]]
  rownames(LM) <- NULL
  LM[1] <- NULL

  #Plot
  if(isTRUE(plot) & !is.null(write)) {
    par(mfrow = c(1,1))
    tiff(paste0(write, '_fragmentacion.tif'), width = 1178, height = 882)
    raster::plot(as(patches, "Spatial"), col = "red")
    raster::plot(CoreA, col = "#1a9641", add = T)
    axis(1)
    axis(2)
    box()
    grid()
    legend("topleft",legend = c("Core Area", "Edge"), fill = c("#1a9641", "red"), cex = 2)
    dev.off()

    p1 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(Area))) +
      ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                     position = "dodge", na.rm = T) +
      ggplot2::labs(x = "log10 (km²)", y = "Frequency", title = "Size") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=20, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20))+
      ggplot2::theme(text = ggplot2::element_text(size = 20),
            axis.text.x = ggplot2::element_text(hjust = 1))
    p2 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(Perimeter))) +
      ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                     position = "dodge", na.rm = T)+
      ggplot2::labs(x = "log10 (km)", y ="Frequency", title = "Perimeter") +
      ggplot2::theme(plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20))+
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     axis.text.x = ggplot2::element_text(hjust = 1))
    p3 <- ggplot2::gplot(data, ggplot2::aes(x = ShapeIndex)) +
      ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                     position = "dodge", na.rm = T) +
      ggplot2::labs(x = "Shape Index", y = "Frequency", title = "Shape Index") +
    theme(plot.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20))+
      ggplot2::theme(text = ggplot2::element_text(size = 20),
            axis.text.x = ggplot2::element_text(hjust = 1))
    p4 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(CA))) +
      ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                     position = "dodge", na.rm = T)+
      ggplot2::labs(x = "log10 (km²)", y = "Frequency", title = "Core Area") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20)) +
      ggplot2::theme(text = ggplot2::element_text(size = 20),
            axis.text.x = ggplot2::element_text(hjust = 1))
    ggplot2::ggarrange(p1, p2, p3, p4)
    ggplot2::ggsave(paste0(write, '_fragStats.tif'), device = "tiff", width = 15,
           height = 11, compression = "lzw")

    } else if (isTRUE(plot) & is.null(write)) {
      p1 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(Area))) +
        ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                       position = "dodge", na.rm = T) +
        ggplot2::labs(x = "log10 (km²)", y = "Frequency", title = "Size") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"),
              axis.title.x = ggplot2::element_text(size = 14),
              axis.title.y = ggplot2::element_text(size = 14))
      p2 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(Perimeter))) +
        ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                       position = "dodge", na.rm = T) +
        ggplot2::labs(x = "log10(km)", y ="Frequency", title = "Perimeter") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=14, face="bold"),
              axis.title.x = ggplot2::element_text(size = 14),
              axis.title.y = ggplot2::element_text(size = 14))
      p3 <- ggplot2::ggplot(data, ggplot2::aes(x = ShapeIndex)) +
        ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10), bins = 10,
                       position = "dodge", na.rm = T) +
        ggplot2::labs(x = "Shape Index", y ="Frequency", title = "Shape Index") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"),
              axis.title.x = ggplot2::element_text(size = 14),
              axis.title.y = ggplot2::element_text(size = 14))
      p4 <- ggplot2::ggplot(data, ggplot2::aes(x = log10(CA))) +
        ggplot2::geom_histogram(color = "black", fill = viridis::viridis(10),
                       bins = 10, position = "dodge", na.rm = T) +
        ggplot2::labs(x = "log10 (km²)", y ="Frequency", title = "Core Area")+
        ggplot2::theme(plot.title = ggplot2::element_text(size=14, face ="bold"),
              axis.title.x = ggplot2::element_text(size = 14),
              axis.title.y = ggplot2::element_text(size = 14))

      p5 <- ggplot2::ggarrange(p1, p2, p3, p4)

    }

  ###Write outputs
  if (!is.null(write)) {
    write.csv(LM, paste0(write, '_LandscapeMetrics.csv'))
    write.csv(data, paste0(write, '_PatchMetrics.csv'))
    sf::write_sf(patches, paste0(write, '_PatchMetrics.shp'), delete_layer = TRUE)
  }

  ###Return
  LM <- formattable::formattable(LM,
              align = c("l","c"),
              list(`Indicator Name` = formattable::formatter("span",
                           style = ~ formattable::style(color = "grey", font.weight = "bold"))))

  if (isTRUE(plot)){
    return(list("Summary landscape metrics (Viewer Panel)" = LM,
                "Patch statistics shapefile" = patches,
                "Frequency landscape metric plot" = p5))
  } else {
    return(list("Summary landscape metrics (Viewer Panel)" = LM,
                "Patch statistics shapefile" = patches))
    }
  }

