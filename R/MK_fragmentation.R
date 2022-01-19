#' Fragmentation Statistics
#'
#' Calculate patch and landscape statistics
#' @param patches Object of class sf, sfc, sfg or SpatialPolygons. Individual patches, the shapefile must be in a projected coordinate system.
#' @param edge_distance Numeric. Distance to edge in meters. Default equal 500 m (Haddad et al. 2015)
#' @param min_patch_area Numeric. Minimum patch area. Default equal 100 km2 (Haddad et al. 2015)
#' @param landscape_area Numeric. Total landscape area in km2 (optional). If NULL the total patch area will be used.
#' @param area_unit character. You can set an area unit (e.g., "km2", "cm2", "m2", "ha"; see Makurhini::unit_convert). Default equal to square kilometers "km2".
#' @param perimeter_unit character. You can set a perimeter unit (e.g., "km", "cm", "m", "ha"; see Makurhini::unit_convert). Default equal to kilometers "km".
#' @param plot Logical. Basic histograms and core area - edge map.
#' @param write Character. Write the tables, shapefile and and plots. It's necessary to specify the path and prefix, for example,
#' to save in the path "C:/Folder" with the prefix "Fragmentation": "C:/Folder/Fragmentation".
#' @return
#' Patch and landscape statistics:\cr
#' 1) Patches Area in square kilometers.\cr
#' 2) Number patches.\cr
#' 3) Mean size of patches.\cr
#' 4) Number of patches smaller than 100 km2 or another minimum patch area.\cr
#' 5) Area in square kilometers of patches smaller than 100 km2 or another minimum patch area.\cr
#' 6) Total edge.\cr
#' 7) Edge density.\cr
#' 8) Total core area (km2; considering a distance to edge of 500 m).\cr
#' 9) Total core area.\cr
#' 10) Core percent.\cr
#' 11) Edge percent.\cr
#' 12) Cority index. It is a measure of fragmentation with respect to a distance from the core area of 500 m (Haddad et al., 2015), where a value of 1 indicates a landscape without fragmentation.\cr
#' 13) Shape Index. A simple shape metric that takes values from 1 (perfectly compact) to infinity is derived by dividing the perimeter by the perimeter of a circle of the same area.\cr
#' 14) Fractal dimension.\cr
#' 15) Effective Mesh Size.\cr
#' @references
#' McGarigal, K., S. A. Cushman, M. C. Neel, and E. Ene. 2002. FRAGSTATS: Spatial Pattern Analysis Program for Categorical Maps. Computer software program produced by the authors at the University of Massachusetts, Amherst. Available at the following web site:
#'  \url{www.umass.edu/landeco/research/fragstats/fragstats.html}.\cr
#' Haddad et al. (2015). Science Advances 1(2):e1500052. \url{DOI: 10.1126/sciadv.1500052}.
#' @examples
#' data("vegetation_patches", package = "Makurhini")
#' nrow(vegetation_patches) # Number of patches
#' fragmentation <- MK_Fragmentation(patches = vegetation_patches, edge_distance = 1000, plot = TRUE)
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
MK_Fragmentation <- function(patches, edge_distance = 500, min_patch_area = 100,
                             landscape_area = NULL, area_unit = "km2", perimeter_unit = "km",
                             plot = FALSE, write = NULL){
  if (missing(patches)) {
    stop("error missing shapefile file of patches")
    } else {
    if (is.numeric(patches) | is.character(patches)) {
      stop("error missing shapefile file of patches")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(patches)[1] == "SpatialPolygonsDataFrame") {
    patches <- st_as_sf(patches) %>% st_cast("POLYGON")
  }

  colors <- c("Core" = "#1a9641", "Edge" = "Red")
  vcol <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")

  patches$IdTemp <- 1:nrow(patches)

  ###Patch metrics
  CoreA <- st_buffer(patches, dist = -(edge_distance))
  data <- data.frame(IdTemp = patches$IdTemp,
                     Area = cbind(round(unit_convert(st_area(patches, byid = T), "m2", area_unit), 4)),
                     CA = cbind(round(unit_convert(st_area(CoreA), "m2", area_unit), 4)))
  data$CAPercent <- round((data$CA * 100) / data$Area, 4)
  data$Perimeter <- round(unit_convert(st_length(st_boundary(patches)), "m", perimeter_unit), 3) %>% as.numeric()
  data$EdgePercent <- round((100 - data$CAPercent), 4)
  data$PARA <- round(data$Area / data$Perimeter, 4)
  data$ShapeIndex <- round(( (data$Perimeter / (2 * pi)) * sqrt(data$Area/pi) ), 4)
  data$FRAC <- round((2 * (log(0.25 * data$Perimeter)) ) / log(data$Area), 4)
  patches <- base::merge(patches, data, by = "IdTemp", all = T)
  patches$IdTemp <- NULL

  ###Landscape metrics
  LM <- data.frame(a = round(sum(data$Area, na.rm = TRUE), 4),
                 b = nrow(data),
                 c = round(mean.default(data$Area, na.rm = TRUE), 4),
                 d = length(which(data$Area < min_patch_area)),
                 e = round(sum(data[which(data$Area < min_patch_area),]$Area, na.rm = TRUE) * 100 / sum(data$Area, na.rm = TRUE), 4),
                 f = round(sum(data$Perimeter, na.rm = TRUE), 4),
                 g = round(sum(data$Perimeter, na.rm = TRUE) / sum(data$Area, na.rm = TRUE), 4),
                 h = round(sum(data$CA, na.rm = TRUE), 4),
                 i = round((nrow(data) - length(which(data$CA == 0))) / nrow(data), 4),
                 j = round(mean(data$ShapeIndex, na.rm = TRUE), 4),
                 k = round(mean(data$FRAC, na.rm = TRUE), 4))

  LM_names <- c(paste0("Patch area (", area_unit, ")"), "Number of patches", "Size (mean)",
                "Patches < minimum patch area", "Patches < minimum patch area (%)",
                "Total edge", "Edge density",
                paste0("Total Core Area (", area_unit, ")"), "Cority",
                "Shape Index (mean)",
                "FRAC (mean)", paste0("MESH (", area_unit, ")"))

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

  if(isTRUE("ggplot2" %in% rownames(installed.packages())) &
     isTRUE("ggpubr" %in% rownames(installed.packages()))){
    if(!is.null(write) & isTRUE(plot)) {
      par(mfrow = c(1,1))
      p0 <- ggplot(data = patches) +
        geom_sf(colour = "Red", aes(fill = "Edge"))+
        geom_sf(data = CoreA[which(!st_is_empty(CoreA)),], colour = "#1a9641", aes(fill = "Core"))+
        scale_fill_manual(name = "Legend", values = colors)+
        theme_bw(base_size = 22)

      ggsave(paste0(write, '_CoreEdge.tif'), plot = p0, device = "tiff", width = 15,
             height = 11, compression = "lzw")

      p1 <- ggplot(data, aes(x = log10(data$Area))) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "log10 (km2)", y = "Frequency", title = "Size") +
        theme(plot.title = element_text(size=20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))+
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p2 <- ggplot(data, aes(x = log10(data$Perimeter))) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T)+
        labs(x = "log10 (km)", y ="Frequency", title = "Perimeter") +
        theme(plot.title = element_text(size = 20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))+
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p3 <- ggplot(data, aes(x = data$ShapeIndex)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "Shape Index", y = "Frequency", title = "Shape Index") +
        theme(plot.title = element_text(size = 20, face = "bold"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))+
        theme(text = element_text(size = 20),
              axis.text.x = element_text(hjust = 1))
      p4 <- ggplot(data, aes(x = log10(data$CA))) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T)+
        labs(x = "log10 (km2)", y = "Frequency", title = "Core Area") +
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

      p1 <- ggplot(data, aes(x = log10(data$Area))) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "log10 (km2)", y = "Frequency", title = "Size") +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      p2 <- ggplot(data, aes(x = log10(data$Perimeter))) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "log10(km)", y ="Frequency", title = "Perimeter") +
        theme(plot.title = element_text(size=14, face="bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      p3 <- ggplot(data, aes(x = data$ShapeIndex)) +
        geom_histogram(color = "black", fill = vcol, bins = 10,
                       position = "dodge", na.rm = T) +
        labs(x = "Shape Index", y ="Frequency", title = "Shape Index") +
        theme(plot.title = element_text(size = 14, face = "bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))
      p4 <- ggplot(data, aes(x = log10(data$CA))) +
        geom_histogram(color = "black", fill = vcol,
                       bins = 10, position = "dodge", na.rm = T) +
        labs(x = "log10 (km2)", y ="Frequency", title = "Core Area")+
        theme(plot.title = element_text(size=14, face ="bold"),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14))

      p5 <- suppressWarnings(ggarrange(p1, p2, p3, p4))

    }
  } else {
    message("To make the plots you need to install the packages ggplot2 and ggpubr")
    plot = FALSE
  }

  ###Write outputs
  if (!is.null(write)) {
    write.csv(LM, paste0(write, '_LandscapeMetrics.csv'))
    write.csv(data, paste0(write, '_PatchMetrics.csv'))
    n <- which(names(patches) %in% c("Area", "CA", "CAPercent", "Perimeter", "EdgePercent", "PARA",
                          "ShapeIndex", "FRAC"))
    patches_2 <- patches
    names(patches_2)[n] <- c("Area", "CA", "CA%", "P", "Edge%", "PARA", "ShapeI", "FRAC")
    write_sf(patches_2, paste0(write, '_PatchMetrics.shp'), delete_layer = TRUE)
  }

  ###Return
  LM <- formattable(LM,
              align = c("l","c"),
              list(`Indicator Name` = formatter("span",
                           style = ~ style(color = "grey", font.weight = "bold"))))

  if (isTRUE(plot)){
    base::print(p0)
    base::print(p5)
    return(list("Summary landscape metrics (Viewer Panel)" = LM,
                "Patch statistics shapefile" = patches))
  } else {
    return(list("Summary landscape metrics (Viewer Panel)" = LM,
                "Patch statistics shapefile" = patches))
    }
  }

