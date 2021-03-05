#' Creation of nodefile attribute.
#'
#' Generates a table with the Core ID and the their attribute value (e.g. area in km²), using a shapefile with these values.
#' @param nodes sp, sf object, RasterLaryer or SpatRaster (terra package). It must be in a projected coordinate system.
#' @param id character. If nodes is a shappefile then you must specify the column name with the node ID
#'  in the shapefile data table. If nodes is a raster layer then raster values (Integer) will be taken as "id".
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute in the data table selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit, ?unit_covert
#' compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to "m2".
#' @param restoration character or vector. If nodes is a shappefile then you must specify the name of the column
#' with restoration value. If nodes is a raster layer then must be a numeric vector with restoration values
#' to each node in the raster. Binary values (0,1), where 1 = existing nodes in the landscape, and 0 = a new node
#'  to add to the initial landscape (restored).
#' @param write character. Output folder path if you use the name option, otherwise, place the output path, with the name and extension ".txt"
#' @return nodo file in .txt format
#' @export
#' @importFrom raster res
#' @importFrom sf st_zm st_as_sf
#' @importFrom rgeos gArea
#' @importFrom methods as
#' @importFrom utils write.table
nodesfile <- function(nodes, id, attribute = NULL, area_unit = "m2",
                      restoration = NULL,
                      write = NULL) {
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }
  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(nodes)[1] == "sf" | class(nodes)[1] == "SpatialPolygonsDataFrame"){
    if(missing(id)){
      stop("Error, missing id argument")
    }

    if(class(nodes)[1] == "sf") {
      nodes <- st_zm(nodes)
      nodes <- as(nodes, 'Spatial')
    } else {
      nodes <- st_as_sf(nodes)
      nodes <- st_zm(nodes)
      nodes <- as(nodes, 'Spatial')
    }

    if (is.null(attribute)) {
      nodes$Area <- gArea(nodes, byid = TRUE)
      attribute <- "Area"
      if (area_unit != "m2"){
        nodes$Area <- unit_convert(nodes$Area, "m2", area_unit)
      }
    }

    if (is.null(restoration)) {
      colvalues <- c(which(names(nodes) == id),
                     which(names(nodes) == attribute))
    } else {
      colvalues <- c(which(names(nodes) == id),
                     which(names(nodes) == attribute),
                     which(names(nodes) == restoration))
    }

    nodes <- nodes@data[ ,colvalues]

  } else if (class(nodes)[1] == "RasterLayer" | class(nodes)[1] == "SpatRaster"){
    nres <- unit_convert(res(nodes)[1]^2, "m2", area_unit)
    nodes <- as.data.frame(table(nodes[]))
    nodes$Freq <- nodes$Freq * nres

    if(!is.null(attribute)){
      nodes$Freq <- nodes$Freq * attribute
    }

    if(is.null(restoration)){
      names(nodes) <- c("Id", "attribute")
      colvalues <- 1:2
    } else {
      nodes$restoration <- restoration
      names(nodes) <- c("Id", "attribute", "restoration")
      colvalues <- 1:3
    }

    nodes <- nodes[ , colvalues]

  } else {
    stop("error missing file of nodes. Check the nodes class")
  }

  if(!is.null(levels(nodes$Id))){
    nodes[,1] <- as.numeric(as.character(nodes$Id))
  }

  if(!is.null(write)){
    if (dim(nodes)[1] > 1){
      write.table(nodes, write, sep="\t", row.names = FALSE, col.names = FALSE)
      return(nodes)
    } else {
      stop("Error only one node")
    }
  }

  return(nodes)
}

