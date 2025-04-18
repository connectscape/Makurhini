#'Estimates PC metric
#'
#' @param x list.
#' @param y distance matrix
#' @param p numeric. Probability
#' @param pmedian logical. median (TRUE) or mean(FALSE) dispersal distance
#' @param d numeric. Dispersal distance
#' @param LA numeric. Max. landscape attribute
#' @importFrom igraph graph_from_adjacency_matrix distances E
#' @importFrom sf st_sf st_as_sf st_geometry
#' @keywords internal
get_pc_grid <- function(x, y, p, pmedian = TRUE, d, LA = NULL){
  st_geometry(x) <- NULL; Adj_matr <- y * 0

  #negative kernel density
  if(isFALSE(pmedian)){
    k <- (1 / d); Adj_matr <- exp(-k * y)
  } else {
    Adj_matr <- exp((y * log(p))/d)
  }

  diag(Adj_matr) <- 0; mode(Adj_matr) <- "numeric"

  #adjacency
  graph_nodes <- tryCatch(igraph::graph_from_adjacency_matrix(Adj_matr, mode = "undirected",
                                          weighted = TRUE), error = function(err) err)

  if (inherits(graph_nodes, "error")) {
    stop("graph adjacency error")
  }

  #product of shortest paths
  pij.mat <- tryCatch(igraph::distances(graph_nodes,
                                     weights = -log(E(graph_nodes)$weight)), error = function(err) err)

  if(inherits(pij.mat, "error")) {
    stop("graph shortest.paths error")
  } else {
    pij.mat <- exp(-pij.mat)
  }

  pij.mat[is.infinite(pij.mat)] <- 1000;  PCmat <- outer(x$attribute, x$attribute) * pij.mat
  PCnum <- sum(PCmat)

  if(!is.null(LA)){
    PC <- PCnum / (LA^2); ECPC <- sqrt(PCnum); ECPC_Normalized <- (ECPC * 100)/LA
    DataPC <- data.frame(cbind(Protected.surface = sum(x[,"attribute"]),
                               LA, 'ECA' = ECPC,
                               'ECA.Normalized' = ECPC_Normalized,
                               PC = PC))
  } else {
    ECPC <- sqrt(PCnum); DataPC <- data.frame(cbind(Protected.surface = sum(x[,"attribute"]),
                               LA, 'ECA' = ECPC))
  }
  return(DataPC)
}
