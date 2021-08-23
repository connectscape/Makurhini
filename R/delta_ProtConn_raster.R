#' Delta ProtConn
#'
#' Estimate the contribution of each node to the ProtConn value in the region.
#' @param x object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param y object of class sf, sfc, sfg or SpatialPolygons. The file must have a projected coordinate system.
#' @param distance_thresholds numeric
#' @param probability numeric
#' @param distance list
#' @param resist raster
#' @param works numeric
#' @importFrom igraph graph.adjacency shortest.paths E
#' @importFrom purrr map_df
#' @importFrom raster beginCluster endCluster clusterR reclassify rasterTmpFile values

delta_ProtConn_raster <- function(x = NULL, y=NULL, distance_thresholds = NULL,
                                  probability = 0.5, LA = NULL,
                                  distance = NULL, resist = NULL,
                                  works = NULL){
  x.1 <- x
  x.m <- y
  if(class(x.1)[1] == "sf"){
    x.1$IdTemp <- 1:nrow(x.1)
    x.1 <- x.1[which(x.1$type != "Region"),]

    distance.d <- tryCatch(protconn_dist(x = x.1, id = "IdTemp",
                                         y = distance,
                                         r = NULL,
                                         resistance = resist),
                           error = function(err)err)

    if(inherits(distance.d, "error")){
      stop("error distance. Check topology errors or resistance raster")
    }

    deltasp <- lapply(distance_thresholds, function(d){
      #Matrix probability
      Adj_matr <- distance.d * 0

      #Negative kernel density
      if(is.null(d)){
        d <- mean(distance.d)
      }

      Adj_matr <- exp((distance.d * log(probability))/d)

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

      aiaj <- outer(x.1$attribute, x.1$attribute, FUN = "*")
      aiaj[t,] <- 1
      aiaj[,t] <- 1
      PCmat <- aiaj * pij.mat
      #p1
      PCnum <- sum(PCmat)
      #ECA
      ECA1 <- sqrt(PCnum)
      #ProtConn
      ProtConn1 <- 100 * (ECA1 / LA)
      Prot1 <- 100* sum(x.1$attribute/LA)#No necesariamente coincide con el prot general que tiene un dissolv

      #delta
      delta.1 <- purrr::map_df(1:max(which(x.1$type == "Non-Transboundary")), function(i){
        mat.i <- Adj_matr[-i,-i]
        n.i <- x.1[-i,]
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
        ProtConn.i <- 100 * (ECA.i / LA)
        Prot.i <- sum(attribute.i)/LA *100

        #Deltas
        delta <- data.frame("dProt" = (Prot1 - Prot.i) / Prot1 * 100,
                            "dProtConn" = (ProtConn1 - ProtConn.i) / ProtConn1 * 100,
                            "varProtConn" = ProtConn1 - ProtConn.i)
        return(delta)})
      delta.1$id <- x.1$id2[which(x.1$type == "Non-Transboundary")]

      v <- unique(raster::values(x.m)) %>% as.numeric()
      v <-v[!is.na(v)]
      v <- sort(v)

      m <- matrix(nrow = length(v), ncol=3, byrow=TRUE)

      for(i in 1:length(v)){
          m[i,1] <- if(i ==1){-Inf} else{v[i-1]}
          m[i,2] <- if(i == length(v)){Inf} else {v[i]}
          m[i,3] <- if(v[i] %in% delta.1$id){delta.1$dProtConn[which(delta.1$id == v[i])]} else{NA}
      }

      if(!is.null(works)){
        beginCluster(works)
        delta.2 <- clusterR(x.m, reclassify, args = list(rcl = m, right = FALSE),
                            filename = rasterTmpFile(), overwrite=TRUE)
        endCluster()
      } else {
        delta.2 <- reclassify(x.m, m)
      }
      return(delta.2)
    })
  } else if(class(x.1)[1] == "data.frame"){
    deltasp <- lapply(distance_thresholds, function(d){
      if(!is.null(works)){
        m <- c(-Inf, x.1$id2, 100, x.1$id2, Inf, NA)
        m <- matrix(m, nrow = 2, ncol=3, byrow=TRUE)
        beginCluster(works)
        delta <- clusterR(x.m, reclassify, args = list(rcl = m, right = FALSE),
                            filename = rasterTmpFile(), overwrite=TRUE)
        endCluster()
      } else {
        delta <- reclassify(x.m, c(-Inf, x.1$id2, 100, x.1$id2, Inf, NA))

      }
      return(delta)
    })
  } else {
    deltasp <- lapply(distance_thresholds, function(d){
      delta <- data.frame("dProt" = NA,
                          "dProtConn" = NA,
                          "varProtConn" = NA)
      return(delta)
    })
  }

  return(deltasp)
}
