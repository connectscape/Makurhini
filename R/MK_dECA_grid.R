#' dA and dECA in a regular grid
#'
#' Use the function to compute the Protected Connected (ProtConn), EC, PC or IIC indexes in a regular grid.
#' @param nodes list of objects class sf, SpatialPolygonsDataFrame. Nodes of each time to analyze.
#' The shapefiles must be in a projected coordinate system.
#' @param area_unit character. Attribute area units. You can set an area unit, "Makurhini::unit_covert()" compatible unit
#' ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to hectares "ha".
#' @param region object of class sf, sfc, sfg or SpatialPolygons. Region shapefile, the shapefile must be
#'  in a projected coordinate system.
#' @param grid_param list. Parameters of the grid shapefile, see \link[Makurhini]{get_grid}.Just omit the parameter 'region'. Example,
#' list(grid_pol = NULL, hexagonal = TRUE, grid_id = NULL, cellsize = unit_convert(1000, "km2", "m2"),
#' grid_boundary = FALSE, clip = FALSE, tolerance = NULL).
#' @param distance list. See \link[Makurhini]{distancefile}. Example, list(type= "centroid", resistance = NULL).
#' @param metric character. Choose a connectivity metric: "IIC" considering topologycal distances or "PC"
#' considering maximum product probabilities.
#' @param distance_threshold numeric. Distance threshold to establish connections (crs units, usually meters).
#' @param probability numeric. Probability of direct dispersal between nodes, Default, 0.5,
#'  that is 50 percentage of probability connection. If probability = NULL, then it will be the inverse of the mean dispersal distance
#' for the species (1/α; Hanski and Ovaskainen 2000).
#' @param intern logical. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @param parallel numeric. Specify the number of cores to use for parallel processing, default = NULL. Parallelize the function using furrr package and multiprocess
#'  plan when there are more than ONE transboundary.
#' @references
#' Matt Strimas-Mackey. \url{http://strimas.com/spatial/hexagonal-grids/}.\cr
#' Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the
#' world's ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S. & Torne, J. (2012). Conefor 2.6. Universidad Politecnica de Madrid. Available
#'  at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. (2006). Comparison and development of new graph-based landscape
#'  connectivity indices: towards the priorization of habitat patches and corridors for conservation.
#'  Landscape Ecology, 21(7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity
#' in landscape conservation planning: comparison with existing indices and application to a case study.
#' Landscape and Urban Planning, 83(2-3): 91-103.
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(rgeos)
#' library(sf)
#'
#' # Four times (T1.2, T2.3, T3.4)
#' data("list_forest_patches", package = "Makurhini")
#' data("study_area", package = "Makurhini")
#' class(list_forest_patches)
#'
#' hexagons_dECA <- MK_dECA_grid(nodes = list_forest_patches,
#'                               region = study_area,
#'                               area_unit = "ha",
#'                               metric = "IIC",
#'                               grid_param = list(hexagonal = TRUE,
#'                                                 cellsize = unit_convert(100, "km2", "m2")),
#'                               distance_threshold = 3000,
#'                               probability = 0.5,
#'                               distance = list(type = "centroid"),
#'                               intern = TRUE,
#'                               parallel = NULL)
#' hexagons_dECA
#' plot(hexagons_dECA["T3.4.dECA"], breaks = "quantile")
#' plot(hexagons_dECA["T3.4.Type.Change"], key.pos = 1)
#'}
#' @export
#' @importFrom sf st_as_sf st_area st_intersection
#' @importFrom future plan multiprocess availableCores
#' @importFrom furrr future_map_dfr
#' @importFrom progressr handlers handler_pbcol progressor
#' @importFrom crayon bgWhite white bgCyan
#' @importFrom purrr map_df map_dfc
#' @importFrom rmapshaper ms_explode

MK_dECA_grid <- function(nodes,
                         area_unit = "m2",
                         region = NULL,
                         grid_param = list(grid_pol = NULL, grid_id = NULL, hexagonal = TRUE,
                                           cellsize = NULL, grid_boundary = FALSE,
                                           clip = FALSE, tolerance = NULL),
                         distance = list(type = "centroid"),
                         metric = "IIC",
                         distance_threshold = NULL,
                         probability = NULL,
                         intern = TRUE, parallel = NULL){
  options(warn = -1)

  if(!is.null(parallel)){
    if(!is.numeric(parallel)){
      stop("if you use parallel argument then you need a numeric value")
    }
  }

  message("Step 1. Reviewing parameters")
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (!metric %in% c("IIC", "PC")) {
    stop("Type must be either 'IIC', or 'PC'")
  }

  if (isTRUE(unique(metric == c("IIC", "PC")))) {
    metric = "IIC"
  }

  if (metric == "PC") {
    if (!is.null(probability) & !is.numeric(probability)) {
      stop("error missing probability")
    }
  }

  if (is.null(distance_threshold)) {
    stop("error missing numeric distance threshold(s)")
  }

  message("Step 2. Grid processing")
  base_grid <- get_grid(region = region,
                          grid_pol = grid_param$grid_pol,
                          grid_id = grid_param$grid_id,
                          hexagonal = grid_param$hexagonal,
                          cellsize = grid_param$cellsize,
                          grid_boundary = grid_param$grid_boundary,
                          clip = grid_param$clip,
                          tolerance = grid_param$tolerance)

  if(class(base_grid)[1] != "grid"){
    stop("error making the grid")
  }

  loop <- 1:nrow(base_grid@grid)

  if (isTRUE(intern)) {
    message("Step 3. Processing metric on the grid. Progress estimated:")
    handlers(global = T)
    handlers(handler_pbcol(complete = function(s) crayon::bgYellow(crayon::white(s)),
                           incomplete = function(s) crayon::bgWhite(crayon::black(s)),
                           intrusiveness = 2))
    pb <- progressor(along = loop)
  } else {
    message("Step 3. Processing metric on the grid")
  }

  if (!is.null(parallel)){
    works <- as.numeric(availableCores())-1
    works <- if(parallel > works){works}else{parallel}
    plan(strategy = multiprocess, gc = TRUE, workers = works)

    nodes <- future_map(nodes, function(x){
      if(class(x)[[1]] == "SpatialPolygonsDataFrame"){
        x <- st_as_sf(x)
      }
      return(x)
    })
    result_1 <- tryCatch(future_map_dfr(loop, function(x){
      if (isTRUE(intern)) {
        pb()
      }

      dECA.2 <- data.frame("Time" = paste0(1:(length(nodes)-1), ".", 2:length(nodes)),
                           "dA" = 0,
                           "dECA" = 0,
                           "comparisons" = "NA",
                           "Type.Change" = "NA")

      x.1 <- base_grid@grid[x,]
      LA <- st_area(x.1)
      LA <- unit_convert(LA, "m2", area_unit)

      nodes.1 <- lapply(nodes, function(i){
          i.1 <- suppressWarnings(st_intersection(i, x.1))
          if(nrow(i.1) > 0){
            i.1 <- ms_explode(i.1)
            i.1$IdTemp <- 1:nrow(i.1)
            i.1 <- i.1["IdTemp"]
          } else {
            i.1 <- "NA"
          }
        return(i.1)
      })

      dECA.2 <- map_dfc(1:(nrow(dECA.2)-1), function(i){
         i.1 <- dECA.2[i,]
         names(i.1) <- paste0("T", i.1$Time, ".",names(i.1))
         i.1 <- i.1[,2:ncol(i.1)]

         if(!is.character(nodes.1[[i]])){
           if(is.character(nodes.1[[i+1]])){ #shp-NA
             i.1[,1:2] <- -100
             i.1[,3] <- dECAfun(-100, -100)
             i.1[,4] <- dECAfun2(-100, -100)
           } else {  #shp-shp
             dECA.1 <- tryCatch(MK_dECA(nodes = nodes.1[i:(i+1)], attribute = NULL,
                                        area_unit = area_unit,
                                        distance = distance,
                                        metric = metric,
                                        probability = probability,
                                        distance_thresholds = distance_threshold,
                                        LA = LA,
                                        plot = FALSE, parallel = NULL,
                                        write = NULL, intern = FALSE), error = function(err)err)
             i.1[,1:4] <- dECA.1[2,8:11]
           }

         } else {
           if(is.character(nodes.1[[i+1]])){
             i.1[,1:2] <- NA
           } else {
             i.1[,1:2] <- 100
             i.1[,3] <- dECAfun(100, 100)
             i.1[,4] <- dECAfun2(100, 100)
           }
         }
         return(i.1)
        })

      dECA.2 <- cbind(x.1, dECA.2)
      return(dECA.2)
      }, .progress = intern ), error = function(err)err)

  } else {
    nodes <- lapply(nodes, function(x){
      if(class(x)[[1]] == "SpatialPolygonsDataFrame"){
        x <- st_as_sf(x)
      }
      return(x)
    })
    x=8
    result_1 <- map_dfr(loop, function(x){
      if (isTRUE(intern)) {
        #pb()
        print(x)
      }

      x.1 <- base_grid@grid[x,]

      dECA.1 <- data.frame("Time" = paste0(1:(length(nodes)-1), ".", 2:length(nodes)),
                           "dA" = 0,
                           "dECA" = 0,
                           "comparisons" = "NA",
                           "Type.Change" = "NA")

      LA <- st_area(x.1)
      LA <- unit_convert(LA, "m2", area_unit)
i=nodes[[1]]
i = st_buffer(i, 0)
      nodes.1 <- lapply(nodes, function(i){
        i.1 <- suppressWarnings(st_intersection(i, x.1))
        if(nrow(i.1) > 0){
          i.1 <- ms_explode(i.1)
          i.1$IdTemp <- 1:nrow(i.1)
          i.1 <- i.1["IdTemp"]
        } else {
          i.1 <- "NA"
        }
        return(i.1)
      })

      dECA.2 <- map_dfc(1:(length(nodes.1)-1), function(i){
        i.1 <- dECA.1[i,]
        names(i.1) <- paste0("T", i.1$Time, ".",names(i.1))
        i.1 <- i.1[,2:ncol(i.1)]

        if(!is.character(nodes.1[[i]])){
          if(is.character(nodes.1[[i+1]])){ #shp-NA
            i.1[,1:2] <- -100
            i.1[,3] <- dECAfun(-100, -100)
            i.1[,4] <- dECAfun2(-100, -100)
          } else {  #shp-shp
            dECA.1 <- tryCatch(MK_dECA(nodes = nodes.1[i:(i+1)], attribute = NULL,
                                       area_unit = area_unit,
                                       distance = distance,
                                       metric = metric,
                                       probability = probability,
                                       distance_thresholds = distance_threshold,
                                       LA = LA,
                                       plot = FALSE, parallel = NULL,
                                       write = NULL, intern = FALSE), error = function(err)err)
            i.1[,1:4] <- dECA.1[2,8:11]
          }

        } else {
          if(is.character(nodes.1[[i+1]])){
            i.1[,1:2] <- NA
          } else {
            i.1[,1:2] <- 100
            i.1[,3] <- dECAfun(100, 100)
            i.1[,4] <- dECAfun2(100, 100)
          }
        }
        return(i.1)
      })
      dECA.2 <- cbind(x.1, dECA.2)
      return(dECA.2)
    })
    }

  result_1$IdTemp  <- NULL
  return(result_1)
}
