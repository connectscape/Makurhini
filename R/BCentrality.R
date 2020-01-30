#' Betweenness Centrality metrics
#'
#' Use this function to calculate the BC, BCIIC and BCPC indexes under one or several distance thresholds.
#'
#' @param nodes Object of class sf, sfc, sfg or SpatialPolygons. The shapefile must be in a projected coordinate system.
#' @param id character. Column name with the nodes id. If NULL, then a new temporal id will be generated.
#' @param attribute character. Column name with the nodes attribute. If NULL, then the patch area (ha) will be estimated and used as the attribute.
#' @param distance list. Distance parameters. For example: type, resistance,or tolerance. For "type" choose one of the distances: "centroid" (faster), "edge", "hausdorff-edge",
#' "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have to use the "resistance" argument.
#'   To See more options consult the help function of distancefile().
#' @param metric character. Choose a Betweenness Centrality Metric: "BC" or "BCIIC" considering topologycal distances or "BCPC" considering maximum product probabilities.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "BPC" metric.
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to ha).
#' @param write character. Write output shapefile, example, "C:/ejemplo.shp".
#'
#' @references Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#'  Freeman L.C. 1977. Set of Measures of Centrality Based on Betweenness. Sociometry 40: 35-41.\cr
#'  Bodin, Ö. & Saura, S. 2010. Ranking individual habitat patches as connectivity providers: integrating network analysis and patch removal experiments. Ecological Modelling 221: 2393-2405.
#' @examples
#' ruta <- system.file("extdata", "Habitat_Patches.shp", package = "Makurhini")
#' cores <- sf::read_sf(ruta)#'
#' nrow(cores) #Number of cores
#' #One distance threshold
#' BCentrality(nodes = cores, id = "id", attribute = NULL, distance = list(type = "centroid"), metric = "BCIIC", distance_thresholds = 30000)
#'
#' #Two or more distance thresholds
#' BCentrality(nodes = cores, id = "id", attribute = NULL, distance = list(type = "centroid"), metric = "BCIIC", distance_thresholds = c(10000, 30000))
#' @export
BCentrality <- function(nodes, id, attribute  = NULL,
                        distance = list(type= "centroid", resistance = NULL, tolerance = NULL, ...),
                        metric = c("BC", "BCIIC", "BCPC"), distance_thresholds = NULL,
                        probability = NULL, LA = NULL, write = NULL) {
  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (!metric %in% c("BC", "BCIIC", "BCPC")) {
    stop("Type must be either 'BC', 'BCIIC', or 'BCPC'")
  }

  if (isTRUE(unique(metric == c("BC", "BCIIC", "BCPC")))) {
    metric = "BC"
  }

  if (metric == "BCPC") {
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

  if (class(nodes)[1] == "SpatialPolygonsDataFrame") {
    nodes <- sf::st_as_sf(nodes)
  }

  options(warn = -1)
  ttt.2 <- getwd()
  temp.1 <- paste0(tempdir(), "/TempInputs", sample(1:1000, 1, replace = TRUE))
  dir.create(temp.1, recursive = T)

  if (is.null(id)) {
    nodes$IdTemp <- 1:nrow(nodes)
  } else {
    nodes$IdTemp <- nodes[[id]]
  }

  nodesfile(nodes, id = "IdTemp", attribute = attribute, area_unit = "ha",
            write = paste0(temp.1, "/nodes.txt"))
  distancefile(nodes, id = "IdTemp", type = distance$type,
               tolerance = distance$tolerance, resistance = distance$resistance,
               threshold = distance$threshold, mask = distance$mask,
               distance_unit = distance$distance_unit, distance$geometry_out,
               write = paste0(temp.1, "/Dist.txt"))

  setwd(temp.1)
  if (is.null(distance$threshold)) {
    pairs = "all"
  } else {
    pairs = "notall"
  }

  pb <- dplyr::progress_estimated(length(distance_thresholds), 0)
  BC_metric <- tryCatch(purrr::map(as.list(distance_thresholds), function(x) {
    if (length(distance_thresholds) > 1) {
      pb$tick()$print()
    }
    tab1 <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                       typeconnection = "dist", typepairs = pairs, index = metric,
                       thdist = x, multdist = NULL, conprob = probability,
                       onlyoverall = FALSE, LA = LA, nrestauration = FALSE,
                       prefix = NULL, write = NULL)
    tab1 <- tab1[["node_importances"]]
    if (!is.null(write)) {
      write <- paste0(write, "_d", x, ".shp")
    }

    tab1 <- merge_conefor(data = tab1, pattern = NULL, merge_shape = nodes,
                          id = "IdTemp", write = write, dA = FALSE, var = FALSE)
    tab1$IdTemp <- NULL
    return(tab1) }), error = function(err) err)

  if (inherits(BC_metric, "error")) {
    setwd(ttt.2)
  } else {
    if (length(distance_thresholds) == 1) {
      BC_metric <- BC_metric[[1]]
    } else {
      names(BC_metric) <- paste0("d", distance_thresholds)
    }
    setwd(ttt.2)
  }

  return(BC_metric)
  }
