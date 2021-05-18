#' Make grid
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param hexagonal logical. If TRUE then the geometry grid will be hexagonal othwewise will be square
#' @param cell_area numeric. Grid area in the CRS units of object x.
#' @param clip logical. If TRUE, the new grid will be clipped to the x area.
#' The operation time will be longer the greater the number of vertices in the polygon of the region.
#' If it is a region with many vertices use the argument **tolerance**.
#' @param tolerance numeric. If Reduces the number of vertices in the region polygon.
#' @param grid_boundary logical.If TRUE, the "hexagons" or "squares" in the boundaries will be removed.
#' @importFrom sf st_as_sf st_sf st_cast st_zm st_simplify st_buffer st_intersection st_area st_make_grid st_crop st_transform st_crs
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom raster extent projection raster projection<-
#' @importFrom rmapshaper ms_dissolve ms_clip
#' @export
make_grid <- function(x, hexagonal = TRUE, cell_area,
                      clip = FALSE, tolerance = NULL, grid_boundary = FALSE) {

  . = NULL

  if(isTRUE(clip) & isTRUE(grid_boundary)){
    stop("Select just one, clip or grid_boundary")
  }

  x <- TopoClean(x)

  if (isFALSE(hexagonal)) {
    cell_width <- sqrt(cell_area)

    ext <- as(extent(x) + cell_width, "SpatialPolygons")
    projection(ext) <- projection(x)
    g <- raster(ext, resolution = cell_width)

    g <- as(g, "SpatialPolygons") %>% st_as_sf() %>% st_cast("POLYGON")
    g <- st_transform(g, crs = st_crs(x))
    inters <- over_poly(g, x)
    g2 <- g[which(inters == 1),]

  } else {
    cell_width <- sqrt(2 * cell_area/sqrt(3))
    g <- st_make_grid(x, cellsize = cell_width, square = F)  %>% st_sf() %>% st_cast("POLYGON")
    inters <- over_poly(g, x)
    g2 <- g[which(inters == 1),]
  }


  if(!is.null(tolerance)){
    region_1 <- tryCatch(st_simplify(x, dTolerance = tolerance,
                                     preserveTopology = TRUE), error = function(err) err)

    if (inherits(region_1, "error")) {
      x <- st_buffer(x, dist = 0)
      x <- tryCatch(st_simplify(x, dTolerance = tolerance, preserveTopology = TRUE),
                    error = function(err) err)

      if (inherits(x, "error")) {
        stop("tolerance error, check geometry errors")
      }

    } else {
      x <- region_1
    }
  }


  if (isTRUE(clip)) {
    g3 <- tryCatch(st_crop(g2, x) %>% ms_clip(., clip = x), error = function(err) err)

    if (inherits(g3, "error")) {
      g3 <- tryCatch(st_intersection(g2, y = st_buffer(x, dist = 0)), error = function(err) err)

      if (inherits(g3, "error")) {
        stop("clip error, check geometry errors")
      }

    }
    g3$rmapshaperid <- NULL

  } else if (isTRUE(grid_boundary)){
    a <- round(as.numeric(st_area(g2[1,])), 0)

    g4 <- tryCatch(st_intersection(g2, y = st_buffer(x, dist = 0)), error = function(err) err)
    g4$Area2 <- round(as.numeric(st_area(g4)),0)

    if (inherits(g4, "error")) {
      stop("error MK_selectbyloc")
    }
    g3 <- g4[which(g4$Area2 >= a),]

  } else {
      g3 <- g2
    }

  row.names(g3) <- as.character(1:nrow(g3))
  g3$OBJECTID <- 1:nrow(g3)
  g3 <- g3[,"OBJECTID"]
  return(g3)
}
