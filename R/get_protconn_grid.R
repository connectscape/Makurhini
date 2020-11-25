#'Estimates ProtConn metric and fractions for a grid
#'
#' @param x list.
#' @param y distance matrix
#' @param p numeric. Probability
#' @param pmedian logical. median (TRUE) or mean(FALSE) dispersal distance
#' @param d numeric. Dispersal distance
#' @param LA numeric. Max. landscape attribute
#' @param bound logical. If TRUE then ProtConn bound will be estimated
#' @importFrom igraph graph.adjacency shortest.paths E
#' @importFrom sf st_sf st_as_sf st_geometry
get_protconn_grid <- function(x, y, p, pmedian = TRUE, d, LA = NULL, bound = FALSE){
  x.1 <- x[[1]][which(x[[1]]$type != "Region"),]
  st_geometry(x.1) <- NULL
  area1 <- x.1[["attribute"]][which(x.1[["type"]] == "Non-Transboundary")]
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
  Adj_matr.1 <- Adj_matr[which(row.names(y) %in% as.character(x.1[,1])),
                         which(colnames(y) %in% as.character(x.1[,1]))]

  graph_nodes <- tryCatch(graph.adjacency(Adj_matr.1, mode = "undirected", weighted = TRUE), error = function(err) err)

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
  t <- which(x.1$type == "Transboundary")
  x.1$attribute[t] <- 1
  PCmat <- outer(x.1$attribute, x.1$attribute) * pij.mat

  #p1
  PCnum <- sum(PCmat)

  #P2
  if(!is.null(LA)){
    PC <- PCnum / (LA^2)
    PC_EC <- sqrt(PCnum)
    DataProtconn <- data.frame(cbind(ECA = PC_EC, PC = PC, LA,
                                     Protected.surface = sum(x.1[,"attribute"])))
  } else {
    stop("error LA is NULL")
  }

  #ProtConn Indicators
  DataProtconn$Prot <- 100 * (DataProtconn$Protected.surface / LA)
  DataProtconn$ProtConn <- 100 * (DataProtconn$ECA / LA)
  DataProtconn$ProtUnconn <- DataProtconn$Prot - DataProtconn$ProtConn
  DataProtconn$RelConn <- 100 * (DataProtconn$ProtConn / DataProtconn$Prot)
  DataProtconn$Unprotected <- 100 - DataProtconn$Prot

  #ProtConn fractions
  DataProtconn$ProtConn_Prot <- (((sqrt(sum(area1^2)) / LA) * 100) / DataProtconn$ProtConn) * 100

  if(length(t) == 0){
    ECAt <- DataProtconn$ECA
    DataProtconn$ProtConn_Trans <- 0
  } else {
    x.2 <- x.1[-t,]

    Adj_matr.2 <- Adj_matr.1[-t,-t]

    #adjacency
    graph_nodes <- tryCatch(graph.adjacency(Adj_matr.2, mode = "undirected", weighted = TRUE), error = function(err) err)

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
    PCmat <- outer(x.2$attribute, x.2$attribute) * pij.mat

    #p1
    PCnum <- sum(PCmat)

    #P2
    ECAt <- sqrt(PCnum)
    DataProtconn$ProtConn_Trans <- 100 * ((100 * ((DataProtconn$ECA - ECAt) / LA)) / DataProtconn$ProtConn)
  }

  DataProtconn$ProtConn_Unprot <- 100 - DataProtconn$ProtConn_Prot - DataProtconn$ProtConn_Trans
  if(DataProtconn$ProtConn_Unprot < 0){
    DataProtconn$ProtConn_Unprot <- 0
  }

  #ProtConn_Prot fractions
  ########ProtConn[Prot] fractions
  within1 <- ((sqrt(sum(x[[2]]^2))) / LA) * 100
  within2 <- sqrt(DataProtconn$Protected.surface / sum(x[[2]]))
  within3 <- within1 * within2

  DataProtconn$ProtConn_Within <- 100 * (within3/DataProtconn$ProtConn)
  DataProtconn$ProtConn_Contig <- 100 - DataProtconn$ProtConn_Within
  if(DataProtconn$ProtConn_Contig < 0){
    DataProtconn$ProtConn_Contig <- 0
  }

  #Land
  DataProtconn$ProtConn_Within_land <- DataProtconn$ProtConn_Within * (DataProtconn$ProtConn / 100)
  DataProtconn$ProtConn_Contig_land <- DataProtconn$ProtConn_Contig * (DataProtconn$ProtConn / 100)
  DataProtconn$ProtConn_Unprot_land <- DataProtconn$ProtConn_Unprot * (DataProtconn$ProtConn / 100)
  DataProtconn$ProtConn_Trans_land <- DataProtconn$ProtConn_Trans * (DataProtconn$ProtConn / 100)

  if(isTRUE(bound)){
    x.3 <- which(x[[1]]$type == "Region")
    if(length(x.2)>0){
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

      x.3 <- x[[1]]$attribute
      t <- which(x[[1]]$type != "Non-Transboundary")
      x.3[t] <- 1
      PCmat <- outer(x.3, x.3) * pij.mat

      PCnum <- sum(PCmat)
      ECAb <- sqrt(PCnum)

      DataProtconn$ProtUnconn_Design <- (100*(ECAb/LA)) - (100*(ECAt/LA))
      DataProtconn$ProtConn_Bound <- DataProtconn$Prot - DataProtconn$ProtUnconn_Design
    } else {
      DataProtconn$ProtUnconn_Design <- NA
      DataProtconn$ProtConn_Bound <- NA
    }

    DataProtconn <- DataProtconn[,c(1:5, 9, 6:7, 19:20, 8, 10:18)]
  }

  return(DataProtconn)
}


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
