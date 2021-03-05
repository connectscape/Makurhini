#'Estimates PC metric
#'
#' @param x list.
#' @param y distance matrix
#' @param p numeric. Probability
#' @param pmedian logical. median (TRUE) or mean(FALSE) dispersal distance
#' @param d numeric. Dispersal distance
#' @param LA numeric. Max. landscape attribute
#' @importFrom igraph graph.adjacency shortest.paths E
#' @importFrom sf st_sf st_as_sf st_geometry
get_pc_grid <- function(x, y, p, pmedian = TRUE, d, LA = NULL){
  st_geometry(x) <- NULL
  Adj_matr <- y * 0

  #negative kernel density
  if(isFALSE(pmedian)){
    k = (1 / d)
    Adj_matr <- exp(-k * y)
  } else {
    k = log(p)/d
    Adj_matr <- exp(k * y)
  }

  diag(Adj_matr) <- 0

  #adjacency
  graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected", weighted = TRUE), error = function(err) err)

  if (inherits(graph_nodes, "error")) {
    stop("graph adjacency error")
  }

  #product of shortest paths
  pij.mat <- tryCatch(shortest.paths(graph_nodes, weights = -log(E(graph_nodes)$weight)), error = function(err) err)

  if(inherits(pij.mat, "error")) {
    stop("graph shortest.paths error")
  } else {
    pij.mat <- exp(-pij.mat)
  }
  PCmat <- outer(x$attribute, x$attribute) * pij.mat

  #p1
  PCnum <- sum(PCmat)

  #P2
  if(!is.null(LA)){
    PC <- PCnum / (LA^2)
    ECPC <- sqrt(PCnum)
    ECPC_Normalized <- (ECPC * 100)/LA
    DataPC <- data.frame(cbind(Protected.surface = sum(x[,"attribute"]),
                               LA, 'ECA' = ECPC,
                               'ECA.Normalized' = ECPC_Normalized,
                               PC = PC))
  } else {
    ECPC <- sqrt(PCnum)
    DataPC <- data.frame(cbind(Protected.surface = sum(x[,"attribute"]),
                               LA, 'ECA' = ECPC))

  }
  return(DataPC)
}
