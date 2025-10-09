#' Dissolve polygones sf or SpatialPolygonsDataFrame
#' @param polygone object of class sf, sfc, sfg or SpatialPolygons
#' @param attribute character. Column name with the dissolve attribute
#' @importFrom magrittr %>%
#' @importFrom sf st_as_sf st_cast st_union
#' @importFrom methods as
#' @keywords internal
pol_dissolve<-function(polygone, attribute){
  . = NULL
  if(class(polygone)[1]!="sf"){
    polygone <- st_as_sf(polygone)
    polygone$Atrib <- polygone[[attribute]]
    polygone <- split(polygone, polygone$Atrib) %>%
      lapply(st_union) %>%
      do.call(c,.) %>%
      st_cast() %>%
      as(.,'Spatial')
  } else {
    polygone <- polygone %>% split(.,paste0(attribute)) %>%
      lapply(st_union) %>%
      do.call(base::c, .)
  }
  return(polygone)
}
