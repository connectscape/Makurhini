#'Estimates ProtConn metric and fractions for a grid
#'
#' @param x list.
#' @param y distance matrix
#' @param p numeric. Probability
#' @param pmedian logical. median (TRUE) or mean(FALSE) dispersal distance
#' @param d numeric. Dispersal distance
#' @param LA numeric. Max. landscape attribute
#' @param bound logical. If TRUE then ProtConn bound will be estimated
#' @importFrom igraph graph.adjacency distances E
#' @importFrom sf st_sf st_as_sf st_geometry
#' @keywords internal
get_protconn_grid <- function(x, y, p, pmedian = TRUE, d, LA = NULL, bound = FALSE){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  diag(y) <- 0; mode(y) <- "numeric"
  y[lower.tri(y, diag = TRUE)] <- NA; dist_nodes <- as.data.frame(as.table(y))
  dist_nodes <- na.omit(dist_nodes); dist_nodes[,1] <- as.numeric(as.character(dist_nodes[,1]))
  dist_nodes[,2] <- as.numeric(as.character(dist_nodes[,2])); names(dist_nodes) <- c("From", "To", "Distance")

  if(is.null(p)){
    k <- (1 / d); dist_nodes$Distance <- -log(exp(-k * dist_nodes$Distance))
  } else {
    dist_nodes$Distance <- -log(exp((dist_nodes$Distance * log(p))/d))
  }

  #Remove "region" x.1
  x.1 <- x[[1]][which(x[[1]]$type != "Region"),]; st_geometry(x.1) <- NULL

  #Area non transb.
  area1 <- x.1[["attribute"]][which(x.1[["type"]] == "Non-Transboundary")]
  graph_nodes <- tryCatch(makegraph(dist_nodes, directed = FALSE), error = function(err) err)

  if (inherits(graph_nodes, "error")) {
    stop("graph adjacency error")
  }

  graph_nodes.1 <- graph_nodes
  if(length(which(x[[1]]$type == "Region"))>0){
    graph_nodes.1$data <- graph_nodes.1$data[-which(graph_nodes.1$data$from %!in% as.character(x.1[[1]]-1)|graph_nodes.1$data$to %!in% as.character(x.1[[1]]-1)),]
    graph_nodes.1$dict <- graph_nodes.1$dict[-which(graph_nodes.1$dict$ref %!in% as.character(x.1[[1]])),]
  }
  #dist_nodes.2 <- dist_nodes[which(dist_nodes$From %in% as.character(x.1[[1]])),]; dist_nodes.2 <- dist_nodes.2[which(dist_nodes.2$To %in% as.character(x.1[[1]])),]
  #graph_nodes.1 <- tryCatch(makegraph(dist_nodes.2, directed = FALSE), error = function(err) err)
  # if (inherits(graph_nodes.1, "error")) {
  #   stop("graph adjacency error")
  # }
  toV <- as.numeric(graph_nodes.1$dict$ref)
  pij.mat <- tryCatch(get_distance_matrix(Graph = graph_nodes.1, from = toV, to = toV, algorithm = "phast"),
                      error = function(err) err)
  if(inherits(pij.mat, "error")) {
    stop("graph shortest.paths error")
  } else {
    pij.mat <- exp(-pij.mat); pij.mat[is.infinite(pij.mat)] <- 0
  }
  invisible(gc())
  tr <- which(x.1$type == "Transboundary"); aiaj <- outer(x.1$attribute, x.1$attribute, FUN = "*")
  aiaj[tr,] <- 1; aiaj[,tr] <- 1; PCmat <- aiaj * pij.mat; PCnum <- sum(PCmat)

  #P2
  if(!is.null(LA)){
    DataProtconn <- data.frame(cbind(ECA = sqrt(PCnum), PC = PCnum / (LA^2), LA,
                                     Protected.surface = sum(area1)))
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
  DataProtconn$ProtConn_Prot <- 100 * ((100 * (sqrt(sum(area1^2)) / LA)) / DataProtconn$ProtConn)

  if(length(tr) == 0){
    ECAn <- DataProtconn$ECA; DataProtconn$ProtConn_Trans <- 0
  } else {
    x.2 <- x.1[-tr,]
    graph_nodes.1$data <- graph_nodes.1$data[-which(graph_nodes.1$data$from %in% (tr-1)|graph_nodes.1$data$to %in% (tr-1)),]
    graph_nodes.1$dict <- graph_nodes.1$dict[-which(graph_nodes.1$dict$ref %in% as.character(tr)),]
    toV <- as.numeric(graph_nodes.1$dict$ref)
    pij.mat <- tryCatch(get_distance_matrix(Graph = graph_nodes.1, from = toV, to = toV, algorithm = "phast"),
                        error = function(err) err)
    if(inherits(pij.mat, "error")) {
      stop("graph shortest.paths error")
    } else {
      pij.mat <- exp(-pij.mat)
    }
    PCmat <- outer(x.2$attribute, x.2$attribute) * pij.mat; PCnum <- sum(PCmat); ECAn <- sqrt(PCnum)
    DataProtconn$ProtConn_Trans <- 100 * ((100 * ((DataProtconn$ECA - ECAn) / LA)) / DataProtconn$ProtConn)
  }

  DataProtconn$ProtConn_Unprot <- 100 - DataProtconn$ProtConn_Prot - DataProtconn$ProtConn_Trans

  if(DataProtconn$ProtConn_Unprot < 0){
    DataProtconn$ProtConn_Unprot <- 0
  }

  #ProtConn_Prot fractions
  ########ProtConn[Prot] fractions
  within1 <- ((sqrt(sum(x[[2]]^2))) / LA) * 100; within2 <- sqrt(DataProtconn$Protected.surface / sum(x[[2]]))
  within3 <- within1 * within2; DataProtconn$ProtConn_Within <- 100 * (within3/DataProtconn$ProtConn)

  if(DataProtconn$ProtConn_Within > 100){
    DataProtconn$ProtConn_Within <- 100
  }

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
    if(length(x.3) > 0){
      toV <- as.numeric(graph_nodes$dict$ref)
      pij.mat <- tryCatch(get_distance_matrix(Graph = graph_nodes, from = toV, to = toV, algorithm = "phast"),
                          error = function(err) err)

      if(inherits(pij.mat, "error")) {
        stop("graph shortest.paths error")
      } else {
        pij.mat <- exp(-pij.mat)
      }

      x.3 <- x[[1]]$attribute; tr <- which(x[[1]]$type != "Non-Transboundary")
      aiaj <- outer(x.3, x.3, FUN = "*"); aiaj[tr,] <- 1; aiaj[,tr] <- 1
      PCmat <- aiaj * pij.mat; PCnum <- sum(PCmat); ECAdesign <- sqrt(PCnum)

      DataProtconn$ProtUnconn_Design <- (100 * (ECAdesign/LA)) - (100 * (DataProtconn$ECA/LA))
      DataProtconn$ProtConn_Bound <- DataProtconn$Prot - DataProtconn$ProtUnconn_Design
    } else {
      DataProtconn$ProtUnconn_Design <- NA
      DataProtconn$ProtConn_Bound <- NA
    }

    DataProtconn <- DataProtconn[,c(1:5, 9, 6:7, 19:20, 8, 10:18)]
  } else {
    DataProtconn <- DataProtconn[,c(1:5, 9, 6:8, 10:18)]
  }

  return(DataProtconn)
}
