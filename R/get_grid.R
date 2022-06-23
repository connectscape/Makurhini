#' ProtConn: input files
#'
#' @param region object of class sf, sfc, sfg or SpatialPolygons
#' @param grid_pol object of class sf, sfc, sfg or SpatialPolygons. Grid hexagones or squares. The
#' shapefile must be in a projected coordinate system.
#' @param grid_id character. Column name of the grid ID.
#' @param hexagonal logical. If FALSE will be a regular grid of "square".
#' @param cellsize numeric. Grid area (square kilometers).
#' @param grid_boundary logical.If TRUE, the Incomplete "hexagons" or "squares" in the boundaries of
#' the region will be discarded
#' @param clip logical. If TRUE, the new grid will be clipped to the region area. The operation time
#' will be longer the greater the number of vertices in the polygon of the region,
#' if it is a region with many vertices use the argument "tolerance".
#' @param tolerance numeric. If "clip" is equal to TRUE reduces the number of vertices in the region polygon.
#' @importFrom sf st_sf st_cast st_buffer st_difference st_area st_geometry st_zm
#' @importFrom magrittr %>%
#' @importFrom rmapshaper ms_dissolve ms_simplify ms_clip
#' @importFrom methods setClass new
#' @keywords internal
get_grid <- function(region = NULL, grid_pol = NULL, grid_id = NULL,
                     hexagonal = TRUE,
                     cellsize = NULL, grid_boundary = FALSE,
                     clip = FALSE, tolerance = NULL){
  class_cache <- new.env(parent = emptyenv())
  MK_grid <- setClass("grid", slots = list(grid = "sf"),
                               where = class_cache)
  if (is.null(grid_pol)) {
    if (is.null(cellsize)) {
      stop("error missing cellsize parameter")
    }

    if (!is.logical(hexagonal)) {
      stop("hexagonal parameter have to be logical value, TRUE = hexagonal, and FALSE = square'")
    }

    if (is.null(region)){
      stop("region parameter have to be an Object of class sf, sfc, sfg or SpatialPolygons'")
    }

    x_grid <- make_grid(x = region, hexagonal = hexagonal,
                        cell_area = cellsize, clip = clip,
                        tolerance = tolerance, grid_boundary = grid_boundary)

  } else {
    x_grid <- TopoClean(grid_pol);x_grid$IdTemp <- 1:nrow(x_grid)
  }

  x_grid$IdTemp <- 1:nrow(x_grid);x_grid <- x_grid[, which(names(x_grid) != "geometry")]
  x_grid <- new("grid", grid =  x_grid)
  return(x_grid)
}
