#' Join a connectivity index (dIIC o dPC, and the "intra",  "dflux", "connector") to core areas
#'
#' @param datat Data table, data frame or direction of node_importances.txt files (character).
#' @param pattern character. Prefix (e.g."*AREA.txt" )
#' @param merge_shape object of class sf, sfc, sfg or SpatialPolygons. It has to have the same "id" used to estimate the node importance.txt
#' @param id character. Column name with the core id
#' @param dA logical. If TRUE, the delta attribute and its variance is selected. Default = FALSE.
#' @param var logical. If TRUE the metric and fractions variance is reteined.
#' @param write character. Write the shapefile, example, "C:/ejemplo/sahapefile.shp".
#' @export
#' @importFrom data.table rbindlist
#' @importFrom sf st_as_sf
#' @importFrom sf st_write
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
      data <- rbindlist(lapply(filenames,fread), fill=T)
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

#' Multiple join connectivity index (dIIC o dPC, and the "intra",  "dflux", "connector") to core areas
#' @param filename character. Direction of node_importances.txt files
#' @param merge_shape SpatialPolygons. Core area shapefile. It has to have the same id used to estimate the node importance.txt
#' @param patterns character. Prefixes, e.g., c("^AREA_25KM.*txt","^AREA_50KM.*txt","^AREA_100KM.*txt")
#' @param id character. Column name with the core area id
#' @param foldname character. Folders names, these folders will be created inside the outputdir
#' @param write character. Output direction
#' @param prefix character. An initial prefix used to write the shapefile, e.g., "AREA"
#' @param merge_coneforlist list. A list of four arguments used in merge_conefor() function: list(colnames, dA, var)
merge_conefor_mult<-function(filename,
                       merge_shape,
                       patterns,
                       id,
                       foldname,
                       write, merge_coneforlist,
                       prefix){
  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }
for (i in 1:length(patterns)){
  dir.create(paste(write,foldname[i],"/", sep = ""))
  Saveshapes = paste(paste(write, foldname[i],"/", sep = ""), prefix, foldname[i],".shp", sep = "")
  if (!is.null(merge_conefor)){
    merge_conefor(datat = filename,
                  pattern = patterns[i], merge_shape = merge_shape,
                  id = id, dA = merge_coneforlist$dA,
                  var = merge_coneforlist$var, write = Saveshapes)
  } else {
    merge_conefor(datat = filename,
                  pattern = patterns[i], merge_shape = merge_shape,
                  id = id, var = FALSE, dA = FALSE,
                  write = Saveshapes)
  }
}
  }
