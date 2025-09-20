#' Topological error correction
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param xsimplify logical or numeric.
#' @importFrom sf st_as_sf st_is_valid st_is_empty st_zm
#' @importFrom terra buffer vect
#' @importFrom rmapshaper ms_simplify
#' @keywords internal
TopoClean <- function(x, xsimplify = FALSE) {
  if(class(x)[1] == "SpatialPolygonsDataFrame" | class(x)[1] == "SpatVector") {
    x <- st_as_sf(x)
  }

  names(x)[grep("geom|geometry", names(x))] <- "geometry"; st_geometry(x) <- "geometry"

  if(any(!st_is_valid(x))){
    message("Topology errors detected. Makurhini will try to fix them, but it may take time. You can cancel and repair them in a GIS instead.")
    x$idtemp <- 1:nrow(x)
    Nvalidos <- x[!st_is_valid(x),]
    Nvalidos2 <- st_is_valid(Nvalidos, reason = TRUE) %>% sub("\\[.*", "", .)
    pb <- txtProgressBar(0,length(Nvalidos2), style = 3)
    Nvalidos_list <- lapply(1:nrow(Nvalidos), function(y){
      setTxtProgressBar(pb, y)
      if(grepl("nested", Nvalidos2[y])){
        suppressMessages(sf_use_s2(FALSE))
        pg  <- st_cast(st_boundary(Nvalidos[y,]), "MULTILINESTRING") %>%
          st_polygonize(.)
        poly <- st_collection_extract(pg, "POLYGON") %>%
          st_cast(., "MULTIPOLYGON", warn = FALSE) %>%
          st_make_valid(.)
        suppressMessages(sf_use_s2(TRUE))
        maskp <- poly[order(st_area(poly), decreasing = TRUE),]
        maskp <- maskp[2:nrow(maskp),]
        maskp <- st_buffer(maskp, 0) %>% ms_simplify(., keep = 0.9) %>% st_buffer(., 10)
        poly2 <- ms_erase(poly, maskp) %>% ms_dissolve() %>% ms_explode()
        poly2 <- st_combine(poly2) %>% st_union()
        attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
        poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL
        return(poly2)
      } else if(grepl("outside shell", Nvalidos2[y])){
        suppressMessages(sf_use_s2(FALSE))
        pg  <- st_cast(st_boundary(Nvalidos[y,]), "MULTILINESTRING") %>%
          st_polygonize(.)
        poly <- st_collection_extract(pg, "POLYGON") %>%
          st_cast(., "MULTIPOLYGON", warn = FALSE) %>%
          st_make_valid(.)
        suppressMessages(sf_use_s2(TRUE))

        poly2 <- st_buffer(poly, 0) %>% ms_simplify(., keep = 0.9) %>% st_buffer(., 10)
        poly2 <- st_combine(poly2) %>% st_union()
        attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
        poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL

        if(isFALSE(areas_similares(p1 = Nvalidos[y,], p2 = poly2, tol = 0.10))){
          suppressMessages(sf_use_s2(FALSE))
          pg  <- st_cast(st_boundary(Nvalidos[y,]), "MULTILINESTRING") %>%
            sf::st_node(.) %>% sf::st_polygonize(.)
          x_fix3 <- st_collection_extract(pg, "POLYGON", warn = FALSE) %>%
            st_cast(., "MULTIPOLYGON", warn = FALSE) %>% st_make_valid(.)
          suppressMessages(sf_use_s2(TRUE))

          maskp <- x_fix3[order(st_area(x_fix3), decreasing = TRUE),]
          maskp <- maskp[2:nrow(maskp),]
          maskp <- st_buffer(maskp, 0) %>% ms_simplify(., keep = 0.9) %>% st_buffer(., 10)
          poly2 <- ms_erase(x_fix3, maskp) %>% ms_dissolve() %>% ms_explode()
          poly2 <- st_combine(poly2) %>% st_union()
          attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
          poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL

          if(isFALSE(areas_similares(p1 = Nvalidos[y,], p2 = poly2, tol = 0.10))){
            poly2 <- st_combine(x_fix3) %>% st_union()
            attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
            poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL
            return(poly2)
          } else {
            return(poly2)
          }
        } else {
          return(poly2)
        }
      } else if(grepl("Ring Self-intersection", Nvalidos2[y])){
        suppressMessages(sf_use_s2(FALSE))
        poly <- st_collection_extract(Nvalidos[y,], "POLYGON") %>%
          st_cast(., "MULTIPOLYGON", warn = FALSE) %>%
          st_make_valid(.)
        suppressMessages(sf_use_s2(TRUE))
        poly2 <- ms_dissolve(poly) %>% ms_explode() %>% st_buffer(., 0)
        poly2 <- st_combine(poly2) %>% st_union()
        attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
        poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL

        if(isFALSE(areas_similares(p1 = Nvalidos[y,], p2 = poly2, tol = 0.10))){
          suppressMessages(sf_use_s2(FALSE))
          pg  <- st_cast(st_boundary(Nvalidos[y,]), "MULTILINESTRING") %>%
            sf::st_node(.) %>% sf::st_polygonize(.)
          x_fix3 <- st_collection_extract(pg, "POLYGON", warn = FALSE) %>%
            st_cast(., "MULTIPOLYGON", warn = FALSE) %>% st_make_valid(.)
          suppressMessages(sf_use_s2(TRUE))

          maskp <- x_fix3[order(st_area(x_fix3), decreasing = TRUE),]
          maskp <- maskp[2:nrow(maskp),]
          maskp <- st_buffer(maskp, 0) %>% ms_simplify(., keep = 0.9) %>% st_buffer(., 10)
          poly2 <- ms_erase(x_fix3, maskp) %>% ms_dissolve() %>% ms_explode()
          poly2 <- st_combine(poly2) %>% st_union()
          attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
          poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL

          if(isFALSE(areas_similares(p1 = Nvalidos[y,], p2 = poly2, tol = 0.10))){
            poly2 <- st_combine(x_fix3) %>% st_union()
            attrs <- st_drop_geometry(Nvalidos[y, , drop = FALSE])
            poly2 <- st_as_sf(cbind(attrs, poly2)); poly2$rmapshaperid <- NULL
            return(poly2)
          } else {
            return(poly2)
          }
        } else {
          return(poly2)
        }
      } else {
        poly2 <- st_buffer(Nvalidos[y,], 0)
        return(poly2)
      }
    })
    Nvalidos_list <- purrr::compact(Nvalidos_list) %>% do.call(rbind, .)
    x2 <- x[-Nvalidos$idtemp,]; x <- rbind(x2, Nvalidos_list)
    x <- x[order(x$idtemp, decreasing = FALSE),]; x$idtemp <- NULL
  }

  if(isTRUE(xsimplify) | is.numeric(xsimplify)){
    x <- ms_simplify(x, keep = if(isTRUE(xsimplify)){0.9} else {xsimplify},
                     method = "vis", keep_shapes = TRUE)
  }

  x <- terra::buffer(vect(x), 0) |> st_as_sf(x = _)
  x <- x[which(st_is_valid(x)),]; x <- x[which(!st_is_empty(x)),]
  return(x)
}


#' Compare areas
#'
#' @param p1 Object of class sf, sfc, sfg or SpatialPolygons
#' @param p2 Object of class sf, sfc, sfg or SpatialPolygons
#' @param tol numeric
#' @importFrom sf st_is_longlat st_transform st_area
#' @keywords internal

areas_similares <- function(p1, p2, tol = 0.10) {
  #stopifnot(nrow(p1) == 1, nrow(p2) == 1)
  # Asegúrate de trabajar en un CRS métrico:
  if (st_is_longlat(p1)) p1 <- st_transform(p1, 3857)
  if (st_is_longlat(p2)) p2 <- st_transform(p2, 3857)

  a1 <- abs(sum(as.numeric(st_area(p1))))
  a2 <- abs(sum(as.numeric(st_area(p2))))

  # diferencia relativa respecto al mayor (robusto)
  rel_diff <- abs(a1 - a2) / max(a1, a2)
  return(rel_diff <= tol)
}
