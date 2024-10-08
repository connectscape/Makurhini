#' Get centroid clumb raster
#' @param x raster with unique ids
#' @param parallel numeric. Specify the number of cores to use for parallel processing
#' @param centroid_geometry logical. Return geometry (TRUE) or data.frame(FALSE)
#' @importFrom raster values crs xyFromCell
#' @importFrom purrr map_dfr
#' @importFrom future plan multicore multisession availableCores
#' @importFrom furrr future_map
#' @importFrom sf st_as_sf
#' @keywords internal
rast_clump_points <- function(x, parallel = NULL, centroid_geometry = FALSE){
  x.0 <- raster::values(x); x.0 <- unique(x.0); x.0 <- x.0[which(!is.na(x.0))]
  crsx <- crs(x)

  if(is.null(parallel)){
    x.1 <- map_dfr(x.0, function(y){
      y.1 <- colMeans(xyFromCell(x, which(x[] == y)))
      return(y.1)
    })
  } else {
    numCores <- as.numeric(availableCores())-1; numCores <- if(parallel > numCores){numCores}else{parallel}
    if(.Platform$OS.type == "unix") {
      strat <- future::multicore
    } else {
      strat <- future::multisession
    }
    plan(strategy = strat, gc = TRUE, workers = numCores)
    x.1 <- future_map_dfr(x.0, function(y){
      y.1 <- colMeans(xyFromCell(x, which(x[] == y)))
      return(y.1)
    }, .progress = TRUE)
    close_multiprocess(numCores)
  }

  if(isTRUE(centroid_geometry)){
    x.2 <- tryCatch(st_as_sf(x.1, coords = c("x", "y"),
                             crs = crsx, agr = "constant"), error = function(err) err)
  } else {
    x.2 <- x.1
  }

  return(x.2)
}
