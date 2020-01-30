#' Creation of nodefile to CONEFOR command line
#'
#' Generates a table with the Core ID and the their attribute value (e.g. area in km²), using a shapefile with these values. The table is the nodefile input for CONEFOR.
#' @param nodes sp or sf object. Nodes shapefile (lines, points or polygons).If area_unit is used then the shapefile must be in a projected coordinate system.
#' @param id character. Column name with the node ID in the shapefile data table.
#' @param attribute character. Column name with the attribute in the data table selected for the nodes. If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit, udunits2 package compatible unit (e.g., "km2", "cm2", "ha"). Default equal to square meters "m2".
#' @param restauration character. Name of the column with restauration value. Binary values (0,1), where 1 = existing nodes in the landscape, and 0 = a new node to add to the initial landscape (restored).
#' @param multiple character. Name of the column with the regions names. Use in case of processing nodes of several independent sites at the same time.
#' @param prefix character. Initial prefix, use in case of processing several sites at the same time in CONEFOR command line.
#' @param write character. Output folder path if you use the name option, otherwise, place the output path, with the name and extension ".txt"
#' @return nodo file in .txt format
#' @export
nodesfile <- function(nodes, id, attribute = NULL, area_unit = "m2", restauration = NULL, multiple = NULL,
                         prefix = NULL, write = NULL) {
  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }
  if(missing(id)){
    stop("Error, missing id argument")
  }
  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(nodes)[1] == "sf") {
    nodes <- sf::st_zm(nodes)
    nodes <- as(nodes, 'Spatial')
  } else {
    nodes <- sf::st_as_sf(nodes)
    nodes <- sf::st_zm(nodes)
    nodes <- as(nodes, 'Spatial')
    }

  if (is.null(attribute)) {
    nodes$Area <- rgeos::gArea(nodes, byid = TRUE)
    attribute <- "Area"
    if (area_unit != "m2"){
      nodes$Area <- udunits2::ud.convert(nodes$Area, "m2", area_unit)
    }
  }

  if (is.null(restauration)) {
      colvalues <- c(which(names(nodes) == id),
                     which(names(nodes) == attribute))
    } else {
      colvalues <- c(which(names(nodes) == id),
                     which(names(nodes) == attribute),
                     which(names(nodes) == restauration))
      }

  if (!is.null(multiple)) {
    colvalues <- c(colvalues[1], which(colnames(nodes) == multiple),
                   colvalues[2:length(colvalues)])
    }

  nodes <- nodes@data[ ,colvalues]
  nodes[,1] <- as.integer(nodes[,1])

  if(is.null(multiple)){
    if(is.null(write)){
      return(nodes)
    } else {
      if (dim(nodes)[1] > 1){
        write.table(nodes, write, sep="\t", row.names = FALSE, col.names = FALSE)
        return(nodes)
    } else {
      stop("Error only one node")
    }
      }
    } else {
      x <- unique(nodes[,2])
      for(i in x) {
        multiple_1 <- nodes[nodes[2] == i,]
        if (dim(multiple_1)[1] > 1){
          if (is.null(restauration)) {
            multiple_2 <- c(multiple_1[1], multiple_1[3])
            } else {
              multiple_2 <- c(multiple_1[1], multiple_1[3], multiple_1[4])
            }
          if (!is.null(prefix)) {
            write.table(multiple_2, paste0(write, paste(prefix, i, sep = "_"), ".txt", sep = "."),
                    sep = "\t", row.names = FALSE, col.names = FALSE)
            } else {
              write.table(multiple_2, paste0(write, paste(i, sep = "_"), ".txt", sep = "."),
                    sep = "\t", row.names = FALSE, col.names = FALSE)
              }
          } else {
            stop("Error only one node")
          }
      }
    }
  }

