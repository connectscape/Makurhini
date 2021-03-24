#' Delta ProtConn
#'
#' Estimate the contribution of each node to the ProtConn value in the region.
#' @param x object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param y object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param base_param3 list
#' @import igraph
#' @importFrom purrr map_df
#' @export

delta_ProtConn <- function(x = NULL, y=NULL, base_param3 = NULL){
  x.1 <- x
  x.m <- base_param3[[1]]@nodes[which(base_param3[[1]]@nodes$IdTemp %in% x.1$IdTemp),]

  if(nrow(x.1)>1){
    x.1$IdTemp <- 1:nrow(x.1)
    x.1$type <- "Non-Transboundary"
    x.2 <- y[which(y$type == "Transboundary"),]
    x.2$IdTemp <- 1:nrow(x.2)+ nrow(x.1)
    x.2 <- rbind(x.1[,c("IdTemp", "attribute", "type")], x.2[,c("IdTemp", "attribute", "type")])

    distance.d <- tryCatch(protconn_dist(x = x.2, id = "IdTemp",
                                         y = base_param3[[2]]@distance,
                                         r = NULL,
                                         resistance = base_param3[[4]]),
                           error = function(err)err)
    if(inherits(distance.d, "error")){
      stop("error distance. Check topology errors or resistance raster")
    }

    deltasp <- lapply(base_param3[[2]]@distance_threshold, function(d){
      #Matrix probability
      Adj_matr <- distance.d * 0

      #Negative kernel density
      if(is.null(d)){
        d <- mean(distance.d)
      }

      Adj_matr <- exp((distance.d * log(base_param3[[2]]@probability))/d)

      diag(Adj_matr) <- 0

      mode(Adj_matr) <- "numeric"


      graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected", weighted = TRUE),
                              error = function(err) err)

      if (inherits(graph_nodes, "error")) {
        stop("graph adjacency error delta protconn")
      }

      #product of shortest paths
      pij.mat <- tryCatch(shortest.paths(graph_nodes, weights = -log(E(graph_nodes)$weight)),
                          error = function(err) err)

      if(inherits(pij.mat, "error")) {
        stop("graph shortest.paths error delta protconn")
      } else {
        pij.mat <- exp(-pij.mat)
      }

      t <- which(x.1$type == "Transboundary")

      aiaj <- outer(x.2$attribute, x.2$attribute, FUN = "*")
      aiaj[t,] <- 1
      aiaj[,t] <- 1
      PCmat <- aiaj * pij.mat
      #p1
      PCnum <- sum(PCmat)
      #ECA
      ECA1 <- sqrt(PCnum)
      #ProtConn
      ProtConn1 <- 100 * (ECA1 / base_param3[[5]])
      Prot1 <- 100* sum(x.2$attribute/base_param3[[5]])#No necesariamente coincide con el prot general que tiene un dissolv

      #delta
      delta.1 <- purrr::map_df(1:nrow(x.1), function(i){
        mat.i <- Adj_matr[-i,-i]
        n.i <- x.2[-i,]
        attribute.i <- n.i$attribute
        g.i <- graph.adjacency(mat.i, mode = "undirected",
                               weighted = TRUE)
        mat.i <- shortest.paths(g.i, weights = -log(E(g.i)$weight))
        mat.i <- exp(-mat.i)
        #
        t <- which(n.i$type == "Transboundary")

        aiaj <- outer(attribute.i, attribute.i, FUN = "*")
        aiaj[t,] <- 1
        aiaj[,t] <- 1
        PCmat.i <- aiaj * mat.i

        num.i <- sum(PCmat.i)
        #ECA
        ECA.i <- sqrt(num.i)
        #ProtConn
        ProtConn.i <- 100 * (ECA.i / base_param3[[5]])
        Prot.i <- sum(attribute.i)/base_param3[[5]] *100

        #Deltas
        delta <- data.frame("dProt" = (Prot1 - Prot.i) / Prot1 * 100,
                            "dProtConn" = (ProtConn1 - ProtConn.i) / ProtConn1 * 100,
                            "varProtConn" = ProtConn1 - ProtConn.i)
        return(delta)})
      delta.2 <- cbind(x.m, delta.1)
      delta.2$IdTemp <- NULL
      return(delta.2)
      })
  } else if(nrow(x.1)==1){
    deltasp <- lapply(base_param3[[2]]@distance_threshold, function(d){
      delta <- data.frame("dProt" = 100,
                          "dProtConn" = 100,
                          "varProtConn" = 100)
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
