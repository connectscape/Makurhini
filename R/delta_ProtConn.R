#' Delta ProtConn
#'
#' Estimate the contribution of each node to the ProtConn value in the region.
#' @param x object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param y object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param base_param3 list
#' @importFrom igraph graph_from_adjacency_matrix distances E
#' @importFrom purrr map_df
#' @keywords internal
delta_ProtConn <- function(x = NULL, y=NULL, base_param3 = NULL){
  x.1 <- x
  if(nrow(x.1) > 1){
    x.m <- base_param3[[1]]@nodes[which(base_param3[[1]]@nodes$IdTemp %in% x.1$IdTemp),]
    x.1$IdTemp <- 1:nrow(x.1); x.1$type <- "Non-Transboundary"

    if(length(which(y$type == "Transboundary")) > 0){
      x.2 <- y[which(y$type == "Transboundary"),]; x.2$IdTemp <- 1:nrow(x.2)+ nrow(x.1)
      x.2 <- rbind(x.1[,c("IdTemp", "attribute", "type")], x.2[,c("IdTemp", "attribute", "type")])
    } else {
      x.2 <- x.1[,c("IdTemp", "attribute", "type")]
    }

    dist_nodes <- tryCatch(protconn_dist(x = x.2, id = "IdTemp",
                                         y = base_param3[[2]]@distance,
                                         r = NULL,
                                         resistance = base_param3[[4]]),
                           error = function(err)err)
    if(inherits(dist_nodes, "error")){
      stop("error distance. Check topology errors or resistance raster")
    } else {
      '%!in%' <- function(x,y)!('%in%'(x,y))
      diag(dist_nodes) <- 0; mode(dist_nodes) <- "numeric"
      dist_nodes[lower.tri(dist_nodes, diag = TRUE)] <- NA; dist_nodes <- as.data.frame(as.table(dist_nodes))
      dist_nodes <- na.omit(dist_nodes); dist_nodes[,1] <- as.numeric(as.character(dist_nodes[,1]))
      dist_nodes[,2] <- as.numeric(as.character(dist_nodes[,2])); names(dist_nodes) <- c("From", "To", "Distance")
    }

    deltasp <- lapply(base_param3[[2]]@distance_threshold, function(d){
      if(is.null(base_param3[[2]]@probability)){
        k <- (1 / d); dist_nodes$Distance <- -log(exp(-k * dist_nodes$Distance))
      } else {
        dist_nodes$Distance <- -log(exp((dist_nodes$Distance * log(base_param3[[2]]@probability))/d))
      }

      graph_nodes <- tryCatch(makegraph(dist_nodes, directed = FALSE), error = function(err) err)
      if (inherits(graph_nodes, "error")) {
        stop("graph adjacency error delta protconn")
      }
      toV <- as.numeric(graph_nodes$dict$ref)
      pij.mat <- tryCatch(get_distance_matrix(Graph = graph_nodes, from = toV, to = toV, algorithm = "phast"),
                          error = function(err) err)
      if(inherits(pij.mat, "error")) {
        stop("graph shortest.paths error")
      } else {
        pij.mat <- exp(-pij.mat); pij.mat[is.infinite(pij.mat)] <- 0
      }

      t <- which(x.1$type == "Transboundary"); aiaj <- outer(x.2$attribute, x.2$attribute, FUN = "*")
      aiaj[t,] <- 1; aiaj[,t] <- 1; PCmat <- aiaj * pij.mat
      PCnum <- sum(PCmat); ECA1 <- sqrt(PCnum)
      #ProtConn
      ProtConn1 <- 100 * (ECA1 / base_param3[[5]]); Prot1 <- 100* sum(x.2$attribute/base_param3[[5]])#No necesariamente coincide con el prot general que tiene un dissolv

      #delta
      delta.1 <- purrr::map_df(1:nrow(x.1), function(i){
        graph_nodes.i <- graph_nodes; n.i <- x.2[-i,]; attribute.i <- n.i$attribute
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from == (i-1)|graph_nodes.i$data$to == (i-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref == i),]
        toV <- as.numeric(graph_nodes.i$dict$ref)
        mat.i <- get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast")
        mat.i <- exp(-mat.i); mat.i[is.infinite(mat.i)] <- 0; t <- which(n.i$type == "Transboundary")
        aiaj <- outer(attribute.i, attribute.i, FUN = "*"); aiaj[t,] <- 1; aiaj[,t] <- 1; PCmat.i <- aiaj * mat.i
        num.i <- sum(PCmat.i); ECA.i <- sqrt(num.i)
        #ProtConn
        ProtConn.i <- 100 * (ECA.i / base_param3[[5]]); Prot.i <- sum(attribute.i)/base_param3[[5]] *100
        #Deltas
        delta <- data.frame("dProt" = (Prot1 - Prot.i) / Prot1 * 100,
                            "dProtConn" = (ProtConn1 - ProtConn.i) / ProtConn1 * 100,
                            "varProtConn" = ProtConn1 - ProtConn.i)
        return(delta)})
      delta.2 <- cbind(x.m, delta.1); delta.2$IdTemp <- NULL
      return(delta.2)
    })
  } else if(nrow(x.1)==1){
    deltasp <- lapply(base_param3[[2]]@distance_threshold, function(d){
      delta <- data.frame("dProt" = 100,
                          "dProtConn" = 100,
                          "varProtConn" = 100); delta <- cbind(x.1, delta)
      return(delta)
    })
  } else {
    deltasp <- lapply(base_param3[[2]]@distance_threshold, function(d){
      delta <- data.frame("dProt" = NA,
                          "dProtConn" = NA,
                          "varProtConn" = NA)
      return(delta)
    })
  }

  return(deltasp)
}
