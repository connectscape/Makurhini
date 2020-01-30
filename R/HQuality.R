#' Habitat quality weighted by node area and estandarized

#' @param x RasterLayer. Habitat quality raster, for example the Climate stability index Raster obtained by the function CSI_Index(). Proyected coordinate system equal to y.
#' @param y Nodes of class sf, sfc, sfg or SpatialPolygons. Proyected coordinate system equal to x.
#' @param id Character. Column name with the nodes id.
#' @param area_uni numeric. Area unit to weighte the habitat quality index, udunits2 package compatible unit (e.g., "km2", "cm2", "ha"). If NULL, then it will be equal to square meters "m2".
#' @param by_zone Character. Colname with the zones to which nodes belongs. A shapefile will be exported by zone (e.g., by ANP, by Ecoregion).
#' @param write Character. Provide folder, name and extension of the output shape, e.g., "C:/folder/CSI.shp". If byzone is not NULL you have to provide only the folder output, e.g., "C:/out/".
#' @param SAGA Logical. Optimize the process using SAGA GIS and RSAGA package (see, \url{https://github.com/r-spatial/RSAGA}).
#' @return SpatialPolygons with value of the climatic stability index weighted by the area, standardized and not standardized for all the nodes and by area of interest.
#' @export
#' @import velox
#' @importFrom magrittr %>%
HQuality <- function(x, y, id = NULL, area_unit = NULL,
                  by_zone = NULL, write = NULL, SAGA = TRUE){
  #
  y0 <- y
  if(class(y)[1] == "sf"){
    y = as(y, 'Spatial')
  }

  if(stringr::word(toString(raster::crs(y)),1) == "+proj=longlat") {
    stop("nodes file should be projected, unids in meters")
  }

  if (is.null(id)) {
    y@data$IdTemp <- 1:nrow(y)
  } else {
    y@data$IdTemp <- y@data[[id]]
  }

  if(toString(raster::crs(x)) != toString(raster::crs(y))){
    x <- raster::projectRaster(x, crs = toString(raster::crs(y)))
  }
  #
  x <- velox(x)
  x$crop(ms_dissolve(y))
  x <- x$as.RasterLayer(band = 1)
  #
  y2 <- spex::qm_rasterToPolygons_sp(x)
  y@data$area <- rgeos::gArea(y, byid = T)
  if (!is.null(area_unit)){
    y@data$area <- udunits2::ud.convert(y@data$area, "m2", area_unit) #area
  }
  y <- rgeos::gBuffer(y, byid = TRUE, width = 0)
  #
  if (isFALSE(SAGA)){
    y3 <- sf::st_as_sf(y) %>% sf::st_cast("POLYGON")
    y3$idn <- sample.int(nrow(y3)/2, nrow(y3), replace = T)
    list_y <- split(y3, y3$idn)
    y2 <- sf::st_as_sf(y2) %>% sf::st_cast("POLYGON")


    list_y_err <-tryCatch(purrr::map(list_y, function(x){
      result <- sf::st_intersection(x, y = y2) %>% sf::st_cast("POLYGON")
      return(result)}), error= function(err)err)

    if(inherits(list_y_err, "error")){
      y2 <- sf::st_buffer(y2, dist = 0)
      list_y <- purrr::map(list_y, function(x){
        x <- sf::st_buffer(x, dist = 0)
        result <- sf::st_intersection(x, y = y2) %>% sf::st_cast("POLYGON")
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
  habitat@data$area2 <- rgeos::gArea(habitat, byid=T) #SELECT_CORES
  if (!is.null(area_unit)){
    habitat@data$area2 <- udunits2::ud.convert(habitat@data$area2, "m2", area_unit)
  }
  #
  habitat <- habitat@data #SELECT_ZC
  #
  habitat$quality <- habitat[ ,names(x)]

  ######index
  quality_a <- plyr::ddply(habitat, plyr::.(IdTemp), plyr::summarize, hq1 = Makurhini::mult(area2, quality))
  quality_b <- plyr::ddply(quality_a, plyr::.(IdTemp), plyr::summarize, hq1 = sum(hq1), number = length(IdTemp))
  quality_c <- sp::merge(x = y, y = quality_b, by = "IdTemp")

  quality_a <- plyr::ddply(quality_c@data,.(IdTemp), plyr::summarize, hq2 = Makurhini::div(hq1, area))
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
      quality_c <- sf::st_as_sf(quality_c)
      quality_c <- sf::write_sf(quality_c, write, delete_layer = TRUE)
    } else {
      raster::shapefile(quality_c, write, overwrite = TRUE)
    }
    return(quality_c)

  } else if (!is.null(write) & !is.null(by_zone)){
    a <- which(colnames(quality_c@data) == by_zone)
    colnames(quality_c@data)[a]<-"zone"
    x <- unique(quality_c@data$zone)

    if(class(y0)[1] == "sf"){
      quality_c <- sf::st_as_sf(quality_c)
      for(i in x) {
        sec <- quality_c[quality_c$zone==i,]
        sf::write_sf(sec, paste0(write,"_",i,".shp"), delete_layer = TRUE)
      }
    } else {
      for(i in x) {
        sec <- quality_c[quality_c@data$zone == i,]
        raster::shapefile(sec, paste0(write,"_",i,".shp"), overwrite = TRUE)
      }
    }
  }

  if (class(y0)[1] == "sf"){
    quality_c <- sf::st_as_sf(quality_c)
    return(quality_c)
  } else {
    return(quality_c)
  }
}
