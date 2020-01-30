#' Estimate the integral index of connectivity (IIC) or the probability of connectivity (PC)
#'
#' Use this function to calculate the PC and IIC indexes under one or several distance thresholds.
#' @param nodes Object of class sf, sfc, sfg or SpatialPolygons.The shapefile must be in a projected coordinate system.
#' @param id character. Column name with the nodes id.If NULL, then a new temporal id will be generated.
#' @param attribute character. Column name with the nodes attribute. If NULL, then the patch area (ha) will be estimated and used as the attribute.
#' @param restauration  character. Name of the column with restauration value. Binary values (0,1), where 1 = existing nodes in the landscape, and 0 = a new node to add to the initial landscape (restored).
#' @param distance list. Distance parameters. For example: type, resistance,or tolerance. For "type" choose one of the distances: "centroid" (faster), "edge", "hausdorff-edge",
#' "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have to use the "resistance" argument.
#'   To See more options consult the help function of distancefile().
#' @param metric character. Choose a connectivity metric: "IIC" considering topologycal distances or "PC" considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param overall logical. If TRUE, then the EC index will be added to the result which is transformed into a list. Default equal to FALSE
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to ha).
#' @param write character. Write output shapefile and overall table (if TRUE overall argument).
#'   It is necessary to specify the "Folder direction" + "Initial prefix",  for example, "C:/ejemplo".
#' @references Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. 2006. Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology 21 (7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.
#' @examples ruta <- system.file("extdata", "Habitat_Patches.shp", package = "Makurhini")
#' cores <- sf::read_sf(ruta)
#' nrow(cores)
#' #One distance threshold
#' IIC <- dConnectivity(nodes = cores, id = "id", attribute = NULL, distance = list(type = "centroid"), metric = "IIC", distance_thresholds = 30000)
#' IIC
#' #Two or more distance thresholds
#' PC <- dConnectivity(nodes = cores, id = "id", attribute = NULL, distance = list(type = "centroid"), metric = "PC", probability = 0.5, distance_thresholds = c(5000, 50000))
#' PC
#' @export

dConnectivity <- function(nodes, id = NULL, attribute  = NULL, restauration = NULL,
                        distance = list(type= "centroid", resistance = NULL, tolerance = NULL, ...),
                        metric = c("IIC", "PC"),
                        probability, distance_thresholds = NULL,
                        overall = FALSE, LA = NULL, write = NULL) {
  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (!metric %in% c("IIC", "PC")) {
    stop("Type must be either 'IIC', or 'PC'")
  }

  if (isTRUE(unique(metric == c("IIC", "PC")))) {
    metric = "IIC"
  }

  if (metric == "PC") {
    if (is.null(probability) | !is.numeric(probability)) {
      stop("error missing probability")
    }
  }

  if (is.null(distance_thresholds)) {
    stop("error missing numeric distance threshold(s)")
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(class(nodes)[1] == "sf") {
    nodes <- as(nodes, 'Spatial')
  }
  options(warn = -1)
  ttt.2 <- getwd()
  #
  temp.1 <- paste0(tempdir(), "/TempInputs", sample(1:1000, 1, replace = TRUE))
  dir.create(temp.1, recursive = T)
  #
  if (is.null(id)) {
    nodes@data$IdTemp <- 1:nrow(nodes)
  } else {
    nodes@data$IdTemp <- nodes@data[[id]]
  }

  #
  nodesfile(nodes, id = "IdTemp", attribute = attribute, area_unit = "ha",
            multiple = NULL, restauration = restauration,
            prefix=NULL, write = paste0(temp.1,"/nodes.txt"))

  distancefile(nodes,  id = "IdTemp", type = distance$type, tolerance = distance$tolerance,
               resistance = distance$resistance, threshold = distance$threshold, mask = distance$mask,
               distance_unit = distance$distance_unit, geometry_out = distance$geometry_out,
               write = paste0(temp.1,"/Dist.txt"))

  setwd(temp.1)
  if(!is.null(restauration)){
    rest = TRUE
  } else {
    rest = FALSE
  }
  #
  if (is.null(distance$threshold)) {
    pairs = "all"
  } else {
    pairs = "notall"
  }
  pb <- dplyr::progress_estimated(length(distance_thresholds), 0)

  result <- tryCatch(purrr::map(as.list(distance_thresholds), function(x){
    if(length(distance_thresholds)>1){
      pb$tick()$print()
    }
    tab1 <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                         typeconnection = "dist", typepairs = pairs,
                         index = metric, thdist = x,
                         multdist = NULL, conprob = probability,
                         onlyoverall = FALSE, LA = LA,
                         nrestauration = rest,
                         prefix = NULL, write = NULL)
    dPC <- tab1[["node_importances"]]

    #
    if (!is.null(write)){
      write <- paste0(write, "_d", x,".shp")
    }

    dPC <- merge_conefor(data = dPC, pattern = NULL,
                         merge_shape = nodes, id = "IdTemp",
                         write = write,
                         dA = FALSE, var = FALSE)
    if (is.null(id)) {
      names(dPC)[which(colnames(dPC) == "IdTemp")] <- "id"
      dPC[moveme(names(dPC), "id first")]
    } else {
      dPC$"IdTemp" <- NULL
    }


    if(isTRUE(overall)){
        roverall <- tab1[["overall_indices"]]
        names(roverall) <- c("Index", "Value")
        result_interm <- list()
        result_interm[[1]] <- dPC
        result_interm[[2]] <- roverall
        names(result_interm) <- c(paste0("_node_importances_d",x), paste0("_overall_d", x))
        if (!is.null(write)){
          write.table(roverall, paste0(write, "_overall_d", x), sep="\t",
                      row.names = FALSE, col.names = TRUE)
        }
        remove(dPC, roverall)
      } else {
        result_interm <- dPC
        remove(dPC)
        }
      return(result_interm)}), error = function(err)err)

  if(inherits(result, "error")){
    setwd(ttt.2)
    } else {
      if (isFALSE(overall) && length(distance_thresholds) == 1){
        toDo <- TRUE
        } else {
          toDo <- FALSE
          }

      if (isTRUE(toDo)){
        result <- result[[1]]
        } else {
          names(result) <- paste0("d", distance_thresholds)
          }
      setwd(ttt.2)
      }
  return(result)
}
