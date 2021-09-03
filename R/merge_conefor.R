#' Join a connectivity index (dIIC o dPC, and the "intra",  "dflux", "connector") to core areas
#'
#' @param datat Data table, data frame or direction of node_importances.txt files (character).
#' @param pattern character. Prefix (e.g."*AREA.txt" )
#' @param merge_shape object of class sf, sfc, sfg or SpatialPolygons. It has to have the same "id" used to estimate the node importance.txt
#' @param id character. Column name with the core id
#' @param dA logical. If TRUE, the delta attribute and its variance is selected. Default = FALSE.
#' @param var logical. If TRUE the metric and fractions variance is reteined.
#' @param write character. Write the shapefile, example, "C:/ejemplo/sahapefile.shp".
#' @importFrom data.table rbindlist
#' @importFrom sf st_write st_as_sf
#' @importFrom utils read.csv
#' @keywords internal
merge_conefor <- function(datat = NULL,
                      pattern = NULL, merge_shape = NULL,
                      id = NULL, dA = FALSE, var = FALSE,
                      write = NULL
                      ){
  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if (is.character(datat) | !is.null(pattern)){
      filenames <-list.files(datat, pattern = paste0(pattern), full.names = TRUE)
      data <- rbindlist(lapply(filenames,read.csv), fill=T)
      data <- as.data.frame(data)
      } else {
        data <- as.data.frame(datat)
        }

  #
  if (isTRUE(dA)){
    data[,which(unique(is.na(data)) == TRUE)] <- NULL
    } else {
      data$dA <- NULL
      data$varA <- NULL
      data[,which(unique(is.na(data)) == TRUE)] <- NULL
    }

  if (isFALSE(var)){
    data <- data[which(grepl("var",names(data), fixed = TRUE) == FALSE)]
    }

  merge_shape <- st_as_sf(merge_shape)

  merge_shape <- merge(x = merge_shape, y = data, by.x = paste0(id), by.y = "Node")

    if (!is.null(write)){
      st_write(merge_shape, write, delete_layer = TRUE)
  }
  return(merge_shape)
}
