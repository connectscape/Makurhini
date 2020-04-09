#' Estimate centrality measures.
#'
#' Use this function to calculate centrality measures. under one or several distance thresholds.
#' @param nodes Object of class sf, SpatialPolygonsDataFrame or raster.
#' It must be in a projected coordinate system. If nodes is a raster layer then raster values (Integer) will be taken as "id".
#' @param distance list. Distance parameters. For example: type, resistance,or keep. For "type" choose one
#' of the distances: "centroid" (faster), "edge", "hausdorff-edge",
#' "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have
#' to use the "resistance" argument. To see more arguments see the ?distancefile.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (in meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' If probability = NULL, then it will be the inverse of the mean dispersal distance for the species (1/α; Hanski and Ovaskainen 2000).
#' @param rasterparallel logical. If nodes is "raster" then you can use this argument to assign the metrics values to the nodes raster. It is useful when raster resolution is less than 100 m2.
#' @param write character. Write output shapefile. It is necessary to specify the "Folder direction"
#'  + "Initial prefix",  for example, "C:/ejemplo".
#' @details This function implements Patch-Scale Connectivity or Centrality Measures. Radial measures: degree, strength (using probability argument, for weighted graphs),
#' eigenvector centrality (eigen), and closeness centrality (close). Medial measures: betweenness centrality (BWC),
#' node memberships (cluster), and modularity (modules, using probability argument).
#' The function builds on functions out of Csardi’s ’igraph’ package.
#' @references
#' Borgatti, S. P., & Everett, M. G. (2006). A Graph-theoretic perspective on centrality. Social Networks, 28(4), 466–484. https://doi.org/10.1016/j.socnet.2005.11.005
#' Hanski, I. and Ovaskainen, O. 2000. The metapopulation capacity of a fragmented landscape. Nature 404: 755–758.
#' @export
#' @examples
#' library(Makurhini)
#' library(sf)
#' data("vegetation_patches", package = "Makurhini")
#' nrow(vegetation_patches) # Number of patches
#' #Two distance threshold,
#' centrality_test <- MK_RMCentrality(nodes = vegetation_patches,
#'                                 distance = list(type = "centroid"),
#'                                  distance_thresholds = c(10000, 100000),
#'                                  probability = 0.5,
#'                                  write = NULL)
#'
#' @importFrom dplyr progress_estimated
#' @importFrom magrittr %>%
#' @import sf
#' @importFrom purrr map map_dbl
#' @importFrom igraph graph.adjacency strength evcent closeness betweenness clusters cluster_louvain degree
#' @importFrom raster as.matrix extent raster stack extent<- writeRaster reclassify crs crs<- unique

MK_RMCentrality <- function(nodes,
                            distance = list(type = "centroid"),
                            distance_thresholds = NULL,
                            binary = TRUE,
                            probability = NULL,
                            rasterparallel = FALSE,
                            write = NULL){
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
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

  if(class(nodes)[1] == "SpatialPolygonsDataFrame"| class(nodes)[1] == "sf"){
    if (nrow(nodes) < 2) {
      stop("error, you need more than 2 nodes")
    }
    nodes <- st_as_sf(nodes)
    nodes$IdTemp <- 1:nrow(nodes)
    idT <- "IdTemp"
  } else {
    if(length(raster::unique(nodes)) < 2){
      stop("error, you need more than 2 nodes")
    }
    idT <- NULL
  }

  dist <- distancefile(nodes = nodes,  id = idT, type = distance$type,
                       distance_unit = distance$distance_unit, keep = distance$keep,
                       resistance = distance$resistance,
                       CostFun = distance$CostFun, ngh = distance$ngh,
                       mask = distance$mask,
                       threshold = distance$threshold,
                       geometry_out = distance$geometry_out,
                       bounding_circles = distance$bounding_circles,
                       parallel = distance$parallel,
                       edgeParallel = distance$edgeParallel,
                       pairwise = FALSE,
                       write = NULL)
  if(!is.null(idT)){
    rownames(dist) <- nodes$IdTemp
    colnames(dist) <- nodes$IdTemp
  }

  pb <- progress_estimated(length(distance_thresholds), 0)

  centrality_result <- map(as.list(distance_thresholds), function(x){

    if(length(distance_thresholds) > 1){
      pb$tick()$print()
    }
    nodes.2 <- nodes

    if(isTRUE(binary)){
      Adj_matr <- dist*0
      Adj_matr[dist < x] <- 1
      diag(Adj_matr) <- 0
      graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected"), error = function(err) err)
    } else {
      Adj_matr <- dist*0
      if(is.null(probability)){
        k =(1 / x)
        Adj_matr <- exp(-k * dist)
      } else {
        k = log(probability)/x
        Adj_matr <- exp(k * dist)
      }
      diag(Adj_matr) <- 0
      graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected", weighted = TRUE), error = function(err) err)
    }

    if (inherits(graph_nodes, "error")) {
      stop("graph adjacency error")
    }

    if(isFALSE(binary)){
      metric.strength <- strength(graph_nodes, weights = 1 / E(graph_nodes)$weight)
      metric.eigen <- evcent(graph_nodes, weights = 1 / E(graph_nodes)$weight)
      metric.close <- closeness(graph_nodes, weights = 1 / E(graph_nodes)$weight, normalized = TRUE)
      metric.between <- betweenness(graph_nodes, weights = 1 / E(graph_nodes)$weight)
      metric.Membcomponents <- clusters(graph_nodes)$membership
      metric.modularity <- cluster_louvain(graph_nodes, weights = 1 / E(graph_nodes)$weight)
      modules <- rep(0, nrow(dist))
      for(i in 1:length(metric.modularity)){
        n <- metric.modularity[[i]] %>% as.numeric()
        modules[which(nodes.2$IdTemp %in% n)] <- i
      }

      metric_conn <- cbind(rownames(dist), metric.strength, metric.eigen$vector, metric.close,
                           metric.between, metric.Membcomponents,
                           modules) %>% as.data.frame()
      names(metric_conn) <- c("id", "strength", "eigen", "close", "BWC", "cluster", "modules")

      if(!is.null(idT)){
        nodes.2$strength <- metric.strength
        nodes.2$eigen <- metric.eigen$vector
        nodes.2$close <- metric.close
        nodes.2$BWC <- metric.between
        nodes.2$cluster <- metric.Membcomponents
        nodes.2$modules <- modules

        nodes.2 <- nodes.2[moveme(names(nodes.2), "geometry last")]
        nodes.2$IdTemp <- NULL

        if(!is.null(write)){
          write_sf(nodes.2, paste0(write, "_", "d", x,  ".shp"), delete_layer = TRUE)
        }
      }

    } else {
      metric.degree <- degree(graph_nodes)
      metric.eigen <- evcent(graph_nodes)
      metric.close <- closeness(graph_nodes)
      metric.between <- betweenness(graph_nodes)
      metric.Membcomponents <- clusters(graph_nodes)$membership
      metric.modularity <- cluster_louvain(graph_nodes)
      modules <- rep(0, nrow(dist))

      for(i in 1:length(metric.modularity)){
        n <- metric.modularity[[i]] %>% as.numeric()
        modules[which(nodes.2$IdTemp %in% n)] <- i
      }
      metric_conn <- cbind(rownames(dist), metric.degree, metric.eigen$vector, metric.close,
                           metric.between, metric.Membcomponents,
                           modules) %>% as.data.frame()
      names(metric_conn) <- c("id", "degree", "eigen", "close", "BWC", "cluster", "modules")

      if(!is.null(idT)){
        nodes.2$degree <- metric.degree
        nodes.2$eigen <- metric.eigen$vector
        nodes.2$close <- metric.close
        nodes.2$BWC <- metric.between
        nodes.2$cluster <- metric.Membcomponents
        nodes.2$modules <- modules

        nodes.2 <- nodes.2[moveme(names(nodes.2), "geometry last")]
        nodes.2$IdTemp <- NULL
        if(!is.null(write)){
          write_sf(nodes.2, paste0(write, "_", "d", x,  ".shp"), delete_layer = TRUE)
        }
        }
      }

    if(is.null(idT)){
      rp <- raster::unique(nodes)
      rp <- as.vector(rp)
      rp <- rp[which(!is.na(rp))]

      if(isTRUE(rasterparallel)){
        m <- matrix(nrow = nrow(dist), ncol = 2)
        m[,1] <- rownames(dist) %>% as.numeric()

        plan(strategy = multiprocess)

        r_metric <- future_map(2:ncol(metric_conn), function(c){
          x1 <- metric_conn[, c(1, c)]
          for(i in rp){
            n <- x1[[which(x1[,1]== i), 2]] %>% as.character() %>% as.numeric()
            m[which(m == i),2] <- n
          }
          x1 <- reclassify(nodes, rcl = m)
          return(x1)}, .progress = TRUE)

      } else {
        m <- matrix(nrow = nrow(dist), ncol = 2)
        m[,1] <- rownames(dist) %>% as.numeric()

        r_metric <- map(2:ncol(metric_conn), function(c){
          x1 <- metric_conn[, c(1, c)]
          for(i in rp){
            n <- x1[[which(x1[,1]== i), 2]] %>% as.character() %>% as.numeric()
            m[which(m == i),2] <- n
          }
          x1 <- reclassify(nodes, rcl = m)
          return(x1)})
      }

      nodes.2 <- list()

      for(i in 2:(length(r_metric)+1)){
        nodes.2[[i]] <- r_metric[[i-1]]
      }

      nodes.2[[1]] <- nodes
      nodes.2 <- stack(nodes.2)
      names(nodes.2) <- names(metric_conn)

      if (!is.null(write)){
        n <- names(metric_conn)
        n <- map(as.list(2:length(n)), function(w){
          x1 <- nodes.2[[w]]
          crs(x1) <- crs(nodes.2)
          writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x, ".tif"), overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
        })
      }
    }

    return(nodes.2) })

  if (length(distance_thresholds) == 1){
    centrality_result <- centrality_result[[1]]
  } else {
    names(centrality_result) <- paste0("d", distance_thresholds)
  }

  return(centrality_result)
}

