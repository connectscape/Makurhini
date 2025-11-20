#' Estimate all shortest distance between nodes using weights
#'
#' @param dist_nodes matrix distance or pairwise data.frame distance
#' @param pij_mat matrix probability of settlement
#' @param attr_nodes numeric. Nodes attribute
#' @param metric character. 'PC' or 'IIC'
#' @param probability numeric
#' @param distance_threshold numeric
#' @param igraph_Dijkstra logical. Use or not igraph to implement the Dijkstra Algorithm
#' @param parallel numeric
#' @param loop logical
#' @param min_nodes numeric
#' @param G1 numeric
#' @param pij_min numeric
#' @param return_pij logic
#' @param intern logical
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv object.size
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map_dfc
#' @importFrom igraph distances E as_ids V graph_from_adjacency_matrix
#' @importFrom magrittr %>%
#' @importFrom cppRouting makegraph get_distance_matrix
#' @importFrom stats as.dist
#' @keywords internal
get_sdist <- function(dist_nodes = NULL,
                      pij_mat = NULL,
                      graph_nodes = NULL,
                      attr_nodes = NULL,
                      metric = NULL,
                      probability = NULL,
                      distance_threshold = NULL,
                      igraph_Dijkstra = FALSE,
                      parallel = NULL,
                      loop = TRUE,
                      G1 = 1000,
                      min_nodes = 2000,
                      return_graph = FALSE,
                      pij_min = 0.01,
                      return_pij = TRUE,
                      intern = TRUE){
  #dist_nodes <- dist
  if(isTRUE(igraph_Dijkstra)){
    if(is.null(graph_nodes)){
      Adj_matr <- dist_nodes * 0
      if(metric == "IIC"){
        Adj_matr[dist_nodes < distance_threshold] <- 1
      } else {
        if(is.null(probability)){
          k = (1 / distance_threshold); Adj_matr <- exp(-k * dist_nodes)
        } else {
          Adj_matr <- exp((dist_nodes * log(probability))/distance_threshold)#pij
          if(!is.null(pij_min)){
            Adj_matr[Adj_matr < pij_min] <- 0
          }
        }
      }
      diag(Adj_matr) <- 0; mode(Adj_matr) <- "numeric"; Adj_matr[is.na(dist_nodes)] <- 0

      graph_nodes <- igraph::graph_from_adjacency_matrix(Adj_matr, mode = "undirected",
                                                         weighted = if(metric == "IIC"){NULL}else{TRUE},
                                                         diag = FALSE)
    }
    toV <- as_ids(V(graph_nodes))
    if(length(toV) > min_nodes & isTRUE(loop)){
      seq_n <- seq(1,length(toV), G1)

      if(!is.null(parallel)){
        works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        m <- as.numeric(object.size(graph_nodes))* 0.001

        if(m > 500){
          m <- (m + round(m/3)) *1024^2
          options(future.globals.maxSize= m)
        }

        plan(strategy = strat, gc = TRUE, workers = works)
        smat <- tryCatch(future_map_dfc(seq_n, function(y){
          toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes,
                                                              to = toV.i[which(!is.na(toV.i))],
                                                              weights = if(metric == "IIC"){NULL
                                                              } else {-log(E(graph_nodes)$weight)})
          return(y.1)})|> as.matrix(x = _),
          error = function(err) err); close_multiprocess(works)

        if(inherits(smat, "error")){
          close_multiprocess(works)
          message("error probably due to the lack of memory in the parallel process. Makurhini will try using sequential process or you can change allocated memory (e.g., options(future.globals.maxSize= 2118123520)),
        before run again this function")
          if(length(seq_n) > 1 & isTRUE(intern)){
            pb <- txtProgressBar(0,length(seq_n), style = 3)
          }
          smat <- lapply(seq_n, function(y){
            if(length(seq_n) > 1 & isTRUE(intern)){
              setTxtProgressBar(pb, which(seq_n == y))
            }
            toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes,
                                                                to = toV.i[which(!is.na(toV.i))],
                                                                weights = if(metric == "IIC"){NULL
                                                                } else {-log(E(graph_nodes)$weight)})
            return(y.1)})|> do.call(cbind, args = _)
        }
      } else {
        if(intern){
          pb <- txtProgressBar(0,length(seq_n), style = 3)
        }
        smat <- lapply(seq_n, function(y){
          if(length(seq_n) > 1 & isTRUE(intern)){
            setTxtProgressBar(pb, which(seq_n == y))
          }
          toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes,
                                                              to = toV.i[which(!is.na(toV.i))],
                                                              weights = if(metric == "IIC"){NULL
                                                              } else {-log(E(graph_nodes)$weight)})
          return(y.1)})|> do.call(cbind, args = _)
      }
    } else {
      smat <- tryCatch(igraph::distances(graph_nodes,
                                         weights = if(metric == "IIC"){NULL
                                         } else {-log(E(graph_nodes)$weight)}),
                       error = function(err) err)
    }

    if(isTRUE(return_pij)){
      pij <- Adj_matr
    }
  } else {
    if (is.null(graph_nodes)) {
      if(is.null(pij_mat)){
        if (metric == "IIC") {
          dist_nodes[dist_nodes >= distance_threshold] <- NA; dist_nodes[!is.na(dist_nodes)] <- 1
        } else {
          if(is.null(probability)){
            dist_nodes <- exp(-((1 / distance_threshold)) * dist_nodes)
          } else {
            dist_nodes <- exp((dist_nodes * log(probability))/distance_threshold)#Cost
          }

          if(!is.null(pij_min)){
            dist_nodes[dist_nodes < pij_min] <- 0
          }
          dist_nodes <- -log(dist_nodes)
        }

        if(isTRUE(return_pij)){
          pij <- exp(-dist_nodes)
        }
      } else {
        dist_nodes <- pij_mat
        if(!is.null(pij_min)){
          dist_nodes[dist_nodes < pij_min] <- 0
        }
        dist_nodes <- -log(dist_nodes)
        if(isTRUE(return_pij)){
          pij <- pij_mat; diag(pij) <- 1
        }
      }

      # Convert distance matrix to edge list
      if (is.matrix(dist_nodes)) {
        dist_nodes[lower.tri(dist_nodes, diag = TRUE)] <- NA
        dist_df <- na.omit(as.data.frame(as.table(dist_nodes)))
        names(dist_df) <- c("From", "To", "Distance")
        dist_df$From <- as.numeric(as.character(dist_df$From))
        dist_df$To <- as.numeric(as.character(dist_df$To))
        graph_nodes <- makegraph(dist_df, directed = FALSE)
      } else {
        names(dist_df) <- c("From", "To", "Distance")
        dist_df$From <- as.numeric(as.character(dist_df$From))
        dist_df$To <- as.numeric(as.character(dist_df$To))
        graph_nodes <- makegraph(dist_df, directed = FALSE)
      }
    }

    toV <- as.numeric(graph_nodes$dict$ref); chunk_idx <- split(seq_along(toV), ceiling(seq_along(toV) / G1))

    get_chunk_dist <- function(i) {
      get_distance_matrix(graph_nodes,
                          from = toV,
                          to = toV[chunk_idx[[i]]],
                          algorithm = "phast",
                          allcores = FALSE)
    }

    if (length(toV) >= min_nodes && isTRUE(loop)) {
      if (!is.null(parallel)) {
        works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}
        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }
        m <- as.numeric(object.size(graph_nodes))* 0.001

        if(m > 500){
          m <- (m + round(m/3)) *1024^2
          options(future.globals.maxSize= m)
        }

        plan(strategy = strat, gc = TRUE, workers = works)
        smat <- tryCatch({
          future_map_dfc(seq_along(chunk_idx), get_chunk_dist, .progress = intern)|> as.matrix(x = _)
        }, error = function(e) {
          message("Parallel error. Retrying sequentially...")
          if (intern){pb <- txtProgressBar(min = 0, max = length(chunk_idx), style = 3)}
          lapply(seq_along(chunk_idx), function(i) {
            if (intern) setTxtProgressBar(pb, i)
            get_chunk_dist(i)
          }) |> do.call(cbind, args = _)
        })
        close_multiprocess(works)
      } else {
        if (intern){pb <- txtProgressBar(min = 0, max = length(chunk_idx), style = 3)}
        smat <- lapply(seq_along(chunk_idx), function(i) {
          if (intern) setTxtProgressBar(pb, i)
          get_chunk_dist(i)
        }) |> do.call(cbind, args = _)
      }
    } else {
      smat <- get_distance_matrix(graph_nodes, from = toV, to = toV, algorithm = "phast")
    }

    if (metric == "IIC") {
      smat[smat >= distance_threshold] <- Inf; diag(smat) <- 0
    }
  }

  if(metric == "PC"){
    smat[is.infinite(smat)] <- 1000000000000000000000
    smat <- exp(-smat)

    if(isTRUE(return_pij)){
      pijM <- smat
    }

    if(!is.null(attr_nodes)){
      smat <- outer(attr_nodes, attr_nodes) * smat
    }
  } else {
    smat[is.infinite(smat)] <- 16443701
    if(isTRUE(return_pij)){
        pijM <- smat
    }
    smat <- outer(attr_nodes, attr_nodes) / (1 + smat);
    #mat2.i[which(mat.i == 1e-07)] <- 0
  }

  if(isTRUE(return_graph)){
    if(isTRUE(return_pij)){
      return(list("shortest distances" = smat, "graph" = graph_nodes,
                  "pij" = pij, "pij*" = pijM))
    } else {
      return(list("shortest distances" = smat, "graph" = graph_nodes))
    }
  } else {
    return(smat)
  }
}
