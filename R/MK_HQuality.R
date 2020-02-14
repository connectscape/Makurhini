#' Habitat quality weighted by node area and estandarized

#' @param x RasterLayer. Habitat quality raster. Proyected coordinate system equal to y shapefile.
#' @param y Nodes of class sf, sfc, sfg or SpatialPolygons. Proyected coordinate system equal to x.
#' @param id Character. Column name with the nodes id.
#' @param area_unit numeric. Area unit to weighte the habitat quality index, udunits2 package compatible unit (e.g., "km2", "cm2", "ha"). If NULL, then it will be equal to square meters "m2".
#' @param by_zone Character. Colname with the zones to which nodes belongs. A shapefile will be exported by zone (e.g., by ANP, by Ecoregion).
#' @param write Character. Provide folder, name and extension of the output shape, e.g., "C:/folder/CSI.shp". If byzone is not NULL you have to provide only the folder output, e.g., "C:/out/".
#' @param SAGA Logical. Optimize the process using SAGA GIS and RSAGA package (see, \url{https://github.com/r-spatial/RSAGA}).
#' @return SpatialPolygons with value of the climatic stability index weighted by the area, standardized and not standardized for all the nodes and by area of interest.
#' @export
#' @importFrom velox velox
#' @importFrom methods as
#' @importFrom stringr word
#' @importFrom raster crs projectRaster shapefile
#' @importFrom spex qm_rasterToPolygons_sp
#' @importFrom rgeos gArea gBuffer
#' @import sf
#' @importFrom purrr map
#' @importFrom plyr ddply .
#' @importFrom sp merge
#' @importFrom dplyr summarize
#' @importFrom rmapshaper ms_dissolve
#' @importFrom iterators iter
#' @importFrom foreach foreach %dopar%
MK_HQuality <- function(x, y, id = NULL, area_unit = NULL,
                  by_zone = NULL, write = NULL, SAGA = TRUE){
  #
  y0 <- y
  if(class(y)[1] == "sf"){
    y = as(y, 'Spatial')
  }

  if(word(toString(crs(y)),1) == "+proj=longlat") {
    stop("nodes file should be projected, unids in meters")
  }

  if (is.null(id)) {
    y@data$IdTemp <- 1:nrow(y)
  } else {
    y@data$IdTemp <- y@data[[id]]
  }

  if(toString(crs(x)) != toString(crs(y))){
    x <- projectRaster(x, crs = toString(crs(y)))
  }
  #
  x <- velox(x)
  x$crop(ms_dissolve(y))
  x <- x$as.RasterLayer(band = 1)
  #
  y2 <- qm_rasterToPolygons_sp(x)
  y@data$area <- gArea(y, byid = T)
  if (!is.null(area_unit)){
    y@data$area <- unit_convert(y@data$area, "m2", area_unit) #area
  }
  y <- gBuffer(y, byid = TRUE, width = 0)
  #
  if (isFALSE(SAGA)){
    y3 <- st_as_sf(y) %>% st_cast("POLYGON")
    y3$idn <- sample.int(nrow(y3)/2, nrow(y3), replace = T)
    list_y <- split(y3, y3$idn)
    y2 <- st_as_sf(y2) %>% st_cast("POLYGON")


    list_y_err <- foreach(x = iter(list_y), .errorhandling = 'pass') %dopar%
      {
        result <- st_intersection(x, y = y2) %>% st_cast("POLYGON")
        return(result)
      }


    if(inherits(list_y_err, "error")){
      y2 <- st_buffer(y2, dist = 0)
      list_y <- map(list_y, function(x){
        x <- st_buffer(x, dist = 0)
        result <- st_intersection(x, y = y2) %>% st_cast("POLYGON")
        return(result)})
    } else {
      list_y <- list_y_err
    }

    habitat <- do.call(rbind, list_y)
    habitat$idn <- NULL

  } else {
    dir_tmp <- paste0(tempdir(), "/out")
    if(dir.exists(dir_tmp)){
      unlink(dir_tmp, recursive = TRUE)
    }
    dir.create(dir_tmp, recursive = TRUE)
    save_tmp <- paste0(dir_tmp, "/out.shp")
    habitat <- RSAGA::rsaga.intersect.polygons(layer_a = y, layer_b = y2,
                                               result = save_tmp, load = TRUE)
  }
  #
  habitat@data$area2 <- gArea(habitat, byid=T) #SELECT_CORES
  if (!is.null(area_unit)){
    habitat@data$area2 <- unit_convert(habitat@data$area2, "m2", area_unit)
  }
  #
  habitat <- habitat@data #SELECT_ZC
  #
  habitat$quality <- habitat[ ,names(x)]

  ######index
  quality_a <- ddply(habitat, .(.data$IdTemp), dplyr::summarize, hq1 = mult(.data$area2, .data$quality))
  quality_b <- ddply(quality_a, .(.data$IdTemp), dplyr::summarize, hq1 = sum(.data$hq1), number = length(.data$IdTemp))
  quality_c <- sp::merge(x = y, y = quality_b, by = "IdTemp")

  quality_a <- ddply(quality_c@data, .(.data$IdTemp), dplyr::summarize, hq2 = div(.data$hq1, .data$area))
  quality_c <- sp::merge(x = quality_c, y = quality_a, by = "IdTemp")
  #
  if (is.null(id)) {
    quality_c@data$OBJECTID <- quality_c$IdTemp
    quality_c$IdTemp <- NULL
    quality_c <- quality_c[moveme(names(quality_c), "OBJECTID first")]
  } else {
    quality_c@data$IdTemp <- NULL
  }

  #
  if (!is.null(write) & is.null(by_zone)){
    if(class(y0)[1] == "sf"){
      quality_c <- st_as_sf(quality_c)
      quality_c <- write_sf(quality_c, write, delete_layer = TRUE)
    } else {
      shapefile(quality_c, write, overwrite = TRUE)
    }
    return(quality_c)

  } else if (!is.null(write) & !is.null(by_zone)){
    a <- which(colnames(quality_c@data) == by_zone)
    colnames(quality_c@data)[a]<-"zone"
    x <- unique(quality_c@data$zone)

    if(class(y0)[1] == "sf"){
      quality_c <- st_as_sf(quality_c)
      for(i in x) {
        sec <- quality_c[quality_c$zone==i,]
        write_sf(sec, paste0(write,"_",i,".shp"), delete_layer = TRUE)
      }
    } else {
      for(i in x) {
        sec <- quality_c[quality_c@data$zone == i,]
        shapefile(sec, paste0(write,"_",i,".shp"), overwrite = TRUE)
      }
    }
  }

  if (class(y0)[1] == "sf"){
    quality_c <- st_as_sf(quality_c)
    return(quality_c)
  } else {
    return(quality_c)
  }
}
