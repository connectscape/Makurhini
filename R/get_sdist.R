#' Estimate all shortest distance between nodes using weights
#'
#' @param dist_nodes matrix distance or pairwise data.frame distance
#' @param metric character. 'PC' or 'IIC'
#' @param probability numeric
#' @param distance_threshold numeric
#' @param igraph_Dijkstra logical. Use or not igraph to implement the Dijkstra Algorithm
#' @param parallel numeric
#' @param loop logical
#' @param G1 numeric
#' @param G2 numeric
#' @param intern logical
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map
#' @importFrom igraph graph.adjacency shortest.paths E as_ids V
#' @importFrom magrittr %>%
#' @importFrom cppRouting makegraph get_distance_matrix
#' @importFrom stats as.dist
get_sdist <- function(dist_nodes = NULL,
                      metric = NULL,
                      probability = NULL,
                      distance_threshold = NULL,
                      igraph_Dijkstra = FALSE,
                      parallel = NULL,
                      loop = TRUE,
                      G1 = 1000,
                      G2 = 1000,
                      intern = TRUE){
  if(isTRUE(igraph_Dijkstra)){
    Adj_matr <- dist_nodes * 0

    if(metric == "IIC"){
      Adj_matr[dist_nodes < distance_threshold] <- 1
    } else {
      if(is.null(probability)){
        k = (1 / distance_threshold); Adj_matr <- exp(-k * dist_nodes)
      } else {
        Adj_matr <- exp((dist_nodes * log(probability))/distance_threshold)
      }
    }

    diag(Adj_matr) <- 0; mode(Adj_matr) <- "numeric"
    graph_nodes <- graph.adjacency(Adj_matr, mode = "undirected",
                                            weighted = if(metric == "IIC"){NULL}else{TRUE})

    if(nrow(Adj_matr) > 3000 & isTRUE(loop)){
      toV <- as_ids(V(graph_nodes)); seq_n <- seq(1,length(toV), G1)

      if(!is.null(parallel)){
        works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        plan(strategy = strat, gc = TRUE, workers = works)

        smat <- future_map(seq_n, function(y){
          toV.i <- toV[y:(y + G1-1)];y.1 <- shortest.paths(graph_nodes,
                                                            to = toV.i[which(!is.na(toV.i))],
                                                            weights = if(metric == "IIC"){NULL
                                                            } else {-log(E(graph_nodes)$weight)})
          return(y.1)}); close_multiprocess(works)
      } else {
        if(length(seq_n) > 1 & isTRUE(intern)){
          pb <- txtProgressBar(0,length(seq_n), style = 3)
        }
        smat <- lapply(seq_n, function(y){
          if(length(seq_n) > 1 & isTRUE(intern)){
            setTxtProgressBar(pb, which(seq_n == y))
          }
          toV.i <- toV[y:(y + G1-1)];y.1 <- shortest.paths(graph_nodes,
                                                            to = toV.i[which(!is.na(toV.i))],
                                                            weights = if(metric == "IIC"){NULL
                                                            } else {-log(E(graph_nodes)$weight)})
          return(y.1)})
      }

      smat <- do.call(cbind, smat)

    } else {
      smat <- tryCatch(shortest.paths(graph_nodes,
                                      weights = if(metric == "IIC"){NULL
                                      } else {-log(E(graph_nodes)$weight)}),
                       error = function(err) err)
    }

  } else {
    if(metric == "IIC"){
      correg <- dist_nodes
      if(is.matrix(dist_nodes)){
        dist_nodes[lower.tri(dist_nodes, diag = TRUE)] <- NA; dist_nodes <- as.data.frame(as.table(dist_nodes))
        dist_nodes <- na.omit(dist_nodes); dist_nodes[,1] <- as.numeric(as.character(dist_nodes[,1]))
        dist_nodes[,2]<-as.numeric(as.character(dist_nodes[,2])); names(dist_nodes)<-c("From", "To", "Distance")
      }
      dist_nodes$Distance[which(dist_nodes$Distance < distance_threshold)] <- 1; dist_nodes$Distance[which(dist_nodes$Distance > 1)] <- 0
    } else {
      if(is.matrix(dist_nodes)){
        dist_nodes[lower.tri(dist_nodes, diag = TRUE)] <- NA; dist_nodes <- as.data.frame(as.table(dist_nodes))
        dist_nodes <- na.omit(dist_nodes); dist_nodes[,1] <- as.numeric(as.character(dist_nodes[,1]))
        dist_nodes[,2]<-as.numeric(as.character(dist_nodes[,2])); names(dist_nodes)<-c("From", "To", "Distance")
      }

      if(is.null(probability)){
        k = (1 / distance_threshold); dist_nodes$Distance <- -log(exp(-k * dist_nodes$Distance))
      } else {
        dist_nodes$Distance <- -log(exp((dist_nodes$Distance * log(probability))/distance_threshold))
      }
    }

    toV <- unique(c(dist_nodes$From, dist_nodes$To)); graph_nodes <- makegraph(dist_nodes, directed = FALSE)

    if(length(toV) > 3000 & isTRUE(loop)){
      seq_n <- seq(1,length(toV), G2)
      if(is.null(parallel)){
        if(length(seq_n) > 1 & isTRUE(intern)){
          pb <- txtProgressBar(0,length(seq_n), style = 3)
        }

        smat <- lapply(seq_n, function(y){
          if(length(seq_n) > 1 & isTRUE(intern)){
            setTxtProgressBar(pb, which(seq_n == y))
          }
          toV.i <- toV[y:(y + G2-1)];y.1 <- get_distance_matrix(Graph = graph_nodes,
                                                                  from = toV, to = toV.i[which(!is.na(toV.i))],
                                                                  allcores=FALSE)
          return(y.1)
        })
      } else {
        works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        plan(strategy = strat, gc = TRUE, workers = works)

        smat <- future_map(seq_n, function(y){
          toV.i <- toV[y:(y + G2-1)];y.1 <- get_distance_matrix(Graph = graph_nodes,
                                                                  from = toV, to = toV.i[which(!is.na(toV.i))],
                                                                  allcores=FALSE)
          return(y.1)
        }, .progress = intern); close_multiprocess(works)
      }

      smat <- do.call(cbind, smat)

    } else {
      smat <- get_distance_matrix(Graph = graph_nodes, from = toV, to = toV, allcores = FALSE)
    }

    if(metric == "IIC"){
      correg1 <- smat * 0; correg1[correg < distance_threshold] <- 1
      smat[which(correg1 == 0)] <- Inf; diag(smat) <- 0
    }
  }

  if(metric == "PC"){ smat <- exp(-smat); smat[is.infinite(smat)] <- 0} else {smat[is.infinite(smat)] <- 0.0000001}
  invisible(gc())
  return(smat)
}


