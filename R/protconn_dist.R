#'Estimates nodes distances to ProtConn
#'
#' @param x object of class sf, sfc, sfg or SpatialPolygons.
#' @param id character.
#' @param y list. Distance list.
#' @param r object of class sf, sfc, sfg or SpatialPolygons.
#' @param resistance raster
#' @param resist.units logical. Transform cost units to geographic units by multiplying the
#'  cost by the resolution of the raster
#' @export
#' @importFrom sf st_buffer
#' @importFrom raster crop
#' @importFrom rmapshaper ms_simplify

protconn_dist <- function(x, id, y, r = NULL, resistance = NULL, resist.units = FALSE){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(y$type %in% c("least-cost", "commute-time")){
    if(is.null(resistance)){
      stop("error, you need a resistance raster")
    }
  }


  if(nrow(x[which(x$type == "Non-Transboundary"),])>1){
    distance_base <- tryCatch(distancefile(nodes = x,  id = id, type = y$type,
                                           distance_unit = y$distance_unit, keep = y$keep,
                                           resistance = resistance,
                                           resist.units = y$resist.units,
                                           CostFun = y$CostFun, ngh = y$ngh,
                                           mask = y$mask,
                                           threshold = y$threshold,
                                           geometry_out = y$geometry_out,
                                           bounding_circles = y$bounding_circles,
                                           parallel = y$parallel,
                                           edgeParallel = y$edgeParallel, pairwise = FALSE,
                                           least_cost.java = y$least_cost.java,
                                           cores.java = y$cores.java, ram.java = y$ram.java,
                                           write = NULL), error = function(err)err)

    if (inherits(distance_base, "error")){
      x <- TopoClean(x)
      distance_base <- tryCatch(distancefile(nodes = x,  id = id, type = y$type,
                                             distance_unit = y$distance_unit, keep = y$keep,
                                             resistance = y$resistance,
                                             resist.units = y$resist.units,
                                             CostFun = y$CostFun, ngh = y$ngh,
                                             mask = y$mask,
                                             threshold = y$threshold,
                                             geometry_out = y$geometry_out,
                                             bounding_circles = y$bounding_circles,
                                             parallel = y$parallel,
                                             edgeParallel = y$edgeParallel, pairwise = FALSE,
                                             least_cost.java = y$least_cost.java,
                                             cores.java = y$cores.java, ram.java = y$ram.java,
                                             write = NULL), error = function(err)err)
      if (inherits(distance_base, "error")){
        stop("distance file error")
      }
    }

    #Correction of centroid bound
    bn <- x[[id]][which(x$type == "Region")]

    if(length(bn) > 0){
      if(!is.null(r)){
        if(length(bn) == 1){
          distance_base[which(row.names(distance_base) %in% as.character(bn)),
                        x[[id]][which(x$type != "Non-Transboundary")]] <- 0
          distance_base[x[[id]][which(x$type != "Non-Transboundary")],
                        which(colnames(distance_base) %in% as.character(bn))] <- 0
        } else {
          r <- ms_simplify(r, keep = 0.1,  method = "vis", keep_shapes = TRUE, explode = TRUE)

          for(i in 1:nrow(r)){
            over_nodes <- over_poly(x[which(x$type == "Non-Transboundary"),], y = r[i,], geometry = TRUE)
            distance_base[which(row.names(distance_base) %in% as.character(bn[i])),
                          which(colnames(distance_base) %in% as.character(over_nodes[[id]]))] <- 0

            distance_base[which(row.names(distance_base) %in% as.numeric(over_nodes[[id]])),
                          which(colnames(distance_base) == as.character(bn[i]))] <- 0
          }
        }

      } else {
        stop("error. You need region shapefile (r parameter)")
      }
    }

  } else {
    distance_base <- "No pair nodes"
  }

  return(distance_base)
}
