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
#' @importFrom sf st_as_sf st_buffer st_difference st_area st_geometry st_transform st_crs
#' @importFrom magrittr %>%
#' @importFrom sp disaggregate
#' @importFrom terra vect intersect na.omit aggregate buffer
#' @importFrom rmapshaper ms_dissolve ms_simplify ms_explode
#' @importFrom purrr map_dfr
#' @keywords internal
Protconn_nodes <- function(x, y, buff = NULL, method = "nodes", xsimplify = FALSE,
                           metrunit = "ha", protconn = TRUE, protconn_bound = FALSE,
                           delta = FALSE){
  options(warn = -1)
  . = NULL

  x$id <- 1; x <- x[,"id"]

  if(nrow(y) > 0){
    if(isTRUE(xsimplify)){
      x.0 <- rmapshaper::ms_simplify(x, keep = 0.1,  method = "vis",
                                     keep_shapes = TRUE, explode = TRUE)
      x.1 <- st_buffer(x.0, 0) %>%  ms_dissolve(.)
    } else {
      if(isTRUE(protconn_bound)){
        x.0 <- ms_explode(x)
      }

      x.1 <- x
    }

    y$PROTIDT <- 1:nrow(y)

    '%!in%' <- function(x,y)!('%in%'(x,y))

    y.1 <- over_poly(y, x.1, geometry = TRUE) %>% st_buffer(., 0); y.1$id <- 1

    if(nrow(y.1) > 1){
      f1 <- terra::intersect(vect(y.1[,c("IdTemp", "id",  "geometry")]), vect(x.1[, "geometry"])) %>%
         terra::na.omit(., geom = TRUE)

      f2 <- tryCatch(terra::aggregate(f1, "id") %>% terra::buffer(., 0) %>%
                       as(., 'Spatial') %>%
                       sp::disaggregate(.) %>%
                       st_as_sf(), error = function(err)err)

      if(inherits(f2, "error")){
        f1 <- terra::buffer(f1, 0)
        f2 <- tryCatch(terra::aggregate(f1, "id") %>% terra::buffer(., 1) %>%
                         as(., 'Spatial') %>%
                         sp::disaggregate(.) %>%
                         st_as_sf(), error = function(err)err)
      }

      f2 <- f2[,"geometry"]

      if(isTRUE(protconn)){
        if(nrow(f2) > 1){
          y.1 <- y.1[,"PROTIDT"]; f2$type <- "Non-Transboundary"; f2$attribute <- as.numeric(st_area(f2)) %>%
            unit_convert(., "m2", metrunit)

          #Transboundary
          y.2 <- y[which(y$PROTIDT %!in% y.1$PROTIDT),]

          if(nrow(y.2) > 0){
            f3 <- st_difference(y.1, x.1)

            #Dos metodos
            if(method == "nodes"){
              mask1 <- rmapshaper::ms_simplify(y.1, method = "vis", keep_shapes = TRUE)%>%
                ms_dissolve() %>% st_buffer(., buff)

            } else {
              mask1 <- rmapshaper::ms_simplify(x.1, method = "vis", keep_shapes = TRUE)%>%
                st_buffer(., buff)
            }

            f4 <- over_poly(y.2[,"geometry"], mask1[,"geometry"], geometry = TRUE) %>%
              vect(.) %>% terra::intersect(., vect(mask1[,"geometry"])) %>% st_as_sf(.)

            f5 <- rbind(f3[,"geometry"], f4[,"geometry"]); f5 <- f5[!st_is_empty(f5), ]

            if(nrow(f5) >= 1){
              f5$id <- 1; f5_test <- tryCatch(terra::aggregate(vect(f5), "id") %>% terra::buffer(., 0) %>%
                               as(., 'Spatial') %>%
                               sp::disaggregate(.) %>%
                               st_as_sf(), error = function(err)err)

              if(inherits(f5_test, "error")){
                f5$id <- 1; f5 <- tryCatch(terra::buffer(vect(f5), 0) %>%
                             terra::aggregate(., "id") %>%
                             terra::buffer(., 1) %>%
                             as(., 'Spatial') %>%
                             sp::disaggregate(.) %>%
                             st_as_sf(), error = function(err)err)
              } else {
                f5 <- f5_test
              }

              f5 <- f5[,"geometry"]; f5$attribute <- 0; f5$type <- "Transboundary"

              if(isTRUE(protconn_bound)){
                x.0 <- x.0[,"geometry"]; x.0$type <- "Region"; x.0$attribute <- 0
                f6 <- tryCatch(rbind(f2, f5, x.0), error = function(err)err)

                if(inherits(f6, "error")){
                  if(f6$message == "arguments have different crs"){
                    f2 <- st_transform(f2, st_crs(x.0)); f5 <- st_transform(f5, st_crs(x.0))
                    f6 <- tryCatch(rbind(f2, f5, x.0), error = function(err)err)
                  } else {
                    stop("Review the nodes and nodes input")
                  }
                }
              } else {
                f6 <- rbind(f2, f5)
              }
            } else {
              f6 <- f2
            }
          } else {
            f6 <- f2
          }

          #N1
          f6$OBJECTID<- 1:nrow(f6); f6 <- f6[,c("OBJECTID", "type", "attribute")]

          #N2
          f1b <- tryCatch(as(f1, 'Spatial') %>% sp::disaggregate(.) %>% st_as_sf(.), error = function(err)err)

          if(inherits(f1b, "error")){
            f1b <- map_dfr(1:nrow(f1), function(i){
              i.1 <- as(f1[i,], 'Spatial') %>% sp::disaggregate(.) %>% st_as_sf(.)
              return(i.1)
            })
          }

          f1b$attribute <- as.numeric(st_area(f1b)) %>% unit_convert(., "m2", metrunit)
          f7 <-  f1b["attribute"] %>% st_drop_geometry(.)

          if(isTRUE(delta)){
            f1 <- st_as_sf(f1); f1$attribute <- as.numeric(st_area(f1)) %>% unit_convert(., "m2", metrunit)
            f1$PROTIDT <- NULL; f1$id <- NULL; f8 <- list(nodes_diss = f6, nodes_nondiss = f7, delta = f1)
          } else {
            f8 <- list(nodes_diss = f6, nodes_nondiss = f7)
          }

        } else {
          a1 <- st_area(f2) %>% as.numeric(); a2 <- st_area(x.1) %>% as.numeric()
          if(a1 >= a2){
            f8 <- unit_convert(a2, "m2", metrunit)
          } else {
            f8 <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)
          }
        }
      } else {
        if(nrow(f2) > 1){
          f2$attribute <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)
          f2$OBJECTID<- 1:nrow(f2); f2$type <- "Non-Transboundary"; f8 <- list(f2[,c("OBJECTID", "attribute", "type")])
        } else {
          a1 <- st_area(f2) %>% as.numeric(); a2 <- st_area(x.1) %>% as.numeric()
          if(a1 >= a2){
            f8 <- unit_convert(a2, "m2", metrunit)
          } else {
            f8 <- as.numeric(st_area(f2)) %>% unit_convert(., "m2", metrunit)
          }
        }
      }
    } else if (nrow(y.1) == 1){
      f1 <- terra::intersect(vect(y.1[,"geometry"]), vect(x.1[,"geometry"])) %>%
        terra::na.omit(., geom = TRUE) %>%  st_as_sf(.)

      if(nrow(f1) > 0){
        a1 <- st_area(f1) %>% as.numeric(.); a2 <- st_area(x.1) %>% as.numeric(.)
        if(a1 >= a2){
          f8 <- unit_convert(a2, "m2", metrunit)
        } else {
          f1 <- TopoClean(f1); f8 <- as.numeric(st_area(f1)) %>% unit_convert(., "m2", metrunit)
        }
      } else {
        f8 = "NA"
      }
    } else {
      f8 <- "NA"
    }
  } else {
    f8 <- "NA"
  }
  return(f8)
}
