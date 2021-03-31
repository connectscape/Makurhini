#' Generates ProtConn nodes files
#'
#' @param x object of class sf, sfc, sfg or SpatialPolygons
#' @param y object of class sf, sfc, sfg or SpatialPolygons
#' @param buff numeric
#' @param method character
#' @param xsimplify logical
#' @param metrunit character
#' @param protconn logical. If FALSE then only PC and EC
#' @param protconn_bound logical
#' @param delta logical
#' @export
#' @importFrom sf st_sf st_cast st_buffer st_difference st_area st_geometry
#' @importFrom magrittr %>%
#' @importFrom rmapshaper ms_dissolve ms_simplify ms_clip
Protconn_nodes <- function(x, y, buff = NULL, method = "nodes", xsimplify = FALSE,
                           metrunit = "ha", protconn = TRUE, protconn_bound = FALSE,
                           delta = FALSE){
  options(warn = -1)
  if(nrow(y) > 0){
    if(isTRUE(xsimplify)){
      x.0 <- rmapshaper::ms_simplify(x, keep = 0.1,  method = "vis",
                                     keep_shapes = TRUE, explode = TRUE)
      x.1 <- st_buffer(x.0, 0) %>%  ms_dissolve(.)
    } else {
      x.0 <- st_cast(x, "POLYGON")
      x.1 <- x
    }
    y$PROTIDT <- 1:nrow(y)


    '%!in%' <- function(x,y)!('%in%'(x,y))
    #
    y.1 <- over_poly(y, x.1, geometry = TRUE)
    y.1 <- st_buffer(y.1, 0)

    if(nrow(y.1) > 1){
      f1 <- ms_clip(y.1, x.1)
      f1 <- f1[!st_is_empty(f1), ]
      f2 <- ms_dissolve(st_geometry(f1)) %>% st_buffer(., 0) %>% st_cast("POLYGON") %>% st_sf()

      if(isTRUE(protconn)){
        if(nrow(f2) > 1){
          y.1 <- y.1[,"PROTIDT"]
          f2$type <- "Non-Transboundary"
          f2$attribute <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)

          #Transboundary
          y.2 <- y[which(y$PROTIDT %!in% y.1$PROTIDT),]

          f3 <- st_difference(y.1, x.1)

          #Dos metodos
          if(method == "nodes"){
            mask1 <- rmapshaper::ms_simplify(y.1, method = "vis", keep_shapes = TRUE)%>%
              st_buffer(., buff)

          } else {
            mask1 <- rmapshaper::ms_simplify(x.1, method = "vis", keep_shapes = TRUE)%>%
              st_buffer(., buff)
          }

          f4 <- over_poly(y.2, mask1, geometry = TRUE) %>% ms_clip(., mask1)#
          f5 <- rbind(f3[,"geometry"], f4[,"geometry"])
          f5 <- f5[!st_is_empty(f5), ]

          if(nrow(f5)>=1){
            f5 <- ms_dissolve(st_geometry(f5)) %>% st_buffer(., 0) %>% st_cast("POLYGON") %>% st_sf()
            f5$attribute <- 0
            f5$type <- "Transboundary"

            if(isTRUE(protconn_bound)){
              x.0 <- x.0[,"geometry"]
              x.0$type <- "Region"
              x.0$attribute <- 0
              f6 <- rbind(f2, f5, x.0)
            } else {
              f6 <- rbind(f2, f5)
            }
          } else {
            f6 <- f2
          }

          #N1
          f6$OBJECTID<- 1:nrow(f6)
          f6 <- f6[,c("OBJECTID", "type", "attribute")]

          #N2
          f1$attribute <- as.numeric(st_area(f1)) %>% unit_convert(., "m2", metrunit)
          f7 <- f1[,c("attribute")]
          st_geometry(f7) <- NULL

          if(isTRUE(delta)){
            f1$PROTIDT <- NULL
            f8 <- list(nodes_diss = f6, nodes_nondiss = f7, delta = f1)
          } else {
            f8 <- list(nodes_diss = f6, nodes_nondiss = f7)
          }

        } else {
          f8 <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)
        }

      } else {
        if(nrow(f2) > 1){
          f2$attribute <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)
          f2$OBJECTID<- 1:nrow(f2)
          f2$type <- "Non-Transboundary"
          f8 <- list(f2[,c("OBJECTID", "attribute", "type")])
        } else {
          f8 <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)
        }
      }

    } else if (nrow(y.1) == 1){
      f1 <- ms_clip(y.1, x.1)
      f1 <- f1[!st_is_empty(f1), ]
      f2 <- ms_dissolve(st_geometry(f1)) %>% st_buffer(., 0) %>% st_cast("POLYGON") %>% st_sf()
      f8 <- sum(as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit))
    } else {
      f8 <- "NA"
    }
  } else {
    f8 <- "NA"
  }
  return(f8)
}
