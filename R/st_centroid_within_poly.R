#'st_centroid_within_poly
#'
#'Polygon centroid
#' @param poly object of class sf, sfc, sfg or SpatialPolygons
#' @references
#' \url{https://stackoverflow.com/questions/52522872/r-sf-package-centroid-within-polygon}
#' \url{https://stackoverflow.com/users/3609772/mitch}.
#' @importFrom magrittr %>%
#' @import sf

st_centroid_within_poly <- function(poly){
  centroid <- poly %>% st_centroid()
  in_poly <- st_within(centroid, poly, sparse = F)[[1]]
  if (in_poly) return(centroid)
  centroid_in_poly <- st_point_on_surface(poly)
  return(centroid_in_poly)
  }

#' Number of vertices
#'
#' @param poly Object of class sf, sfc, sfg or SpatialPolygons
#' @import sf
#' @export

num_vert <- function(poly){
  if(class(poly)[1] == "SpatialPolygonsDataFrame") {
    poly <- st_as_sf(poly)
  }
  poly <- st_cast(poly$geometry, "MULTIPOINT")
  poly <- sapply(poly, length)
  poly <- sum(poly)/2
  return(poly)
}

#' over sf
#'
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param y Object of class sf, sfc, sfg or SpatialPolygons
#' @import sf
#' @importFrom magrittr %>%

over_poly <- function(x, y) {
  if(class(x)[1] == "SpatialPolygonsDataFrame") {
    x <- st_as_sf(x)
  }

  if(class(y)[1] == "SpatialPolygonsDataFrame") {
    y <- st_as_sf(y)
  }
  over_result <- sapply(st_intersects(x, st_geometry(y)), function(z) if (length(z)==0) NA_integer_ else z[1])
  return(over_result)
}

#' Make grid
#' @param x Object of class sf, sfc, sfg or SpatialPolygons
#' @param type character. c("square", "hexagonal")
#' @param cell_area numeirc. Grid area (square kilometers).
#' @param clip logical. If TRUE, the new grid will be clipped to the region area. The operation time will be longer the greater the number of vertices in the polygon of the region, if it is a region with many vertices use the argument "tolerance".
#' @param tolerance numeric. If Reduces the number of vertices in the region polygon.
#' @param grid_boundary logical.If TRUE, the Incomplete "hexagons" or "squares" in the boundaries of the region will be discarded
#' @import sf
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom raster extent projection raster
#' @importFrom sp HexPoints2SpatialPolygons spsample
#' @importFrom rmapshaper ms_dissolve
make_grid <- function(x, type, cell_area,
                      clip = FALSE, tolerance = NULL, grid_boundary = NULL) {
  if (class(x)[1] == "sf") {
    x <- st_zm(x)
    x <- as(x, "Spatial")
  }

  if (!type %in% c("square", "hexagonal")) {
    stop("Type must be either 'square' or 'hexagonal'")
  }

  if (type == "square") {
    cell_width <- sqrt(cell_area)
  } else if (type == "hexagonal") {
    cell_width <- sqrt(2 * cell_area/sqrt(3))
  }

  ext <- as(extent(x) + cell_width, "SpatialPolygons")
  projection(ext) <- projection(x)

  if (type == "square") {
    g <- raster(ext, resolution = cell_width)
    g <- as(g, "SpatialPolygons") %>% st_as_sf() %>%
      st_cast("POLYGON")
  } else if (type == "hexagonal") {
    g <- spsample(ext, type = "hexagonal", cellsize = cell_width,
                  offset = c(0, 0))
    g <- HexPoints2SpatialPolygons(g, dx = cell_width) %>%
      st_as_sf() %>% st_cast("POLYGON")
  } else {
    stop("Select a grid_type = c('square', 'hexagonal')")
  }

  x2 <- st_as_sf(x) %>% st_cast("POLYGON")

  if(num_vert(x2) > 100) {
    region_1 <- tryCatch(st_simplify(x2, dTolerance = if(is.null(tolerance)){0}else{tolerance},
                                     preserveTopology = TRUE), error = function(err) err)

    if (inherits(region_1, "error")) {
      x2 <- st_buffer(x2, dist = 0)
      x2 <- st_simplify(x2, dTolerance = if(is.null(tolerance)){0}else{tolerance}, preserveTopology = TRUE)
    } else {
      x2 <- region_1
    }
  }


  x2 <- tryCatch(x2 %>% st_cast("POLYGON") %>% ms_dissolve() %>%
                   st_buffer(dist = 0), error = function(err) err)

  if (inherits(x2, "error")) {
    x2 <- tryCatch(x2 %>% st_buffer(., dist = 0) %>%
                     ms_dissolve() %>% st_buffer(., dist = 0), error = function(err) err)
  }

  if (isTRUE(clip)) {
    if (num_vert(x2) > 100) {
        region_1 <- tryCatch(st_simplify(x2, dTolerance = if(is.null(tolerance)){0}else{tolerance},
          preserveTopology = TRUE), error = function(err) err)

        if (inherits(region_1, "error")) {
          x2 <- st_buffer(x2, dist = 0)
          x2 <- st_simplify(x2, dTolerance = if(is.null(tolerance)){0}else{tolerance}, preserveTopology = TRUE)
        } else {
          x2 <- region_1
        }
      }

    g2 <- tryCatch(st_intersection(g, y = x2) %>% st_cast("POLYGON"),
                   error = function(err) err)

    if (inherits(g2, "error")) {
      g2 <- tryCatch(st_intersection(g, y = st_buffer(x2,dist = 0)) %>%
                       st_cast("POLYGON"), error = function(err) err)
    }

  } else {
    a <- unique(format(round(as.numeric(st_area(g)), 0), scientific = F))

    g2 <- MK_selectbyloc(g, x2)

    if(isTRUE(grid_boundary)){
      if (num_vert(x2) > 100) {
        region_1 <- tryCatch(st_simplify(x2, dTolerance = if(is.null(tolerance)){0}else{tolerance},
                                         preserveTopology = TRUE), error = function(err) err)
        if (inherits(region_1, "error")) {
          x2 <- st_buffer(x2, dist = 0)
          x2 <- st_simplify(x2, dTolerance =if(is.null(tolerance)){0}else{tolerance}, preserveTopology = TRUE)
        } else {
          x2 <- region_1
        }
      }

      g3 <- tryCatch(st_intersection(g2, y = x2) %>% st_cast("POLYGON"),
                     error = function(err) err)

      if (inherits(g3, "error")) {
        g3 <- tryCatch(st_intersection(g2, y = st_buffer(x2,dist = 0)) %>%
                         st_cast("POLYGON"), error = function(err) err)
      }

      g3$area <- format(round(as.numeric(st_area(g3)),0), scientific = F)
      nn <- which(g3$area == a)
      if(length(nn) > 0){
        g2 <- g3[which(g3$area == a), ]
      } else {
        g2 <- g3
      }
    }
  }

  row.names(g2) <- as.character(1:nrow(g2))
  g2$rmapshaperid <- NULL
  g2$area <- NULL
  return(g2)
}
