#' Estimate the integral index of connectivity (IIC) or the probability of connectivity (PC)
#'
#' Use this function to calculate the PC and IIC indexes under one or several distance thresholds.
#' @param nodes Object of class sf, SpatialPolygonsDataFrame, raster or data frame.
#' If spatial layer, then it must be in a projected coordinate system. If nodes is a raster layer then raster values (Integer)
#' will be taken as "id". If nodes is a data frame, then it must have at least two columns, first column with nodes "id" and second with the attribut.
#' If you use the restoration argument, then you must have a dataframe with three columns where the third column is equal the restauration values
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#'   If nodes is a data frame then it must have two columns where second column is the attribute.
#' @param area_unit character. If attribute is NULL you can set an area unit (e.g., "km2", "cm2", "ha";
#'  see Makurhini::unit_convert). Default equal to hectares "ha".
#' @param restauration character or vector. If nodes is a shappefile then you must specify the name of the column
#' with restauration value. If nodes is a raster layer then must be a numeric vector with restauration values
#' to each node in the raster. Binary values (0,1), where 1 = existing nodes in the landscape, and 0 = a new node
#'  to add to the initial landscape (restored).
#'  If nodes is a data frame then it must have three columns where third column has restauration values.
#' @param distance list or distance matrix. If it is a list, then it must be a set of distance parameters. For example: type, resistance,or keep. For "type" choose one
#' of the distances: "centroid" (faster), "edge", "hausdorff-edge",
#' "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have
#' to use the "resistance" argument. To see more arguments see the ?distancefile.
#' If it is a distance matrix, then the number of columns and rows must be equal the number of nodes and "id".
#' @param metric character. Choose a connectivity metric: "IIC" considering topologycal distances or "PC" considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' If probability = NULL, then it will be the inverse of the mean dispersal distance for the species (1/α; Hanski and Ovaskainen 2000).
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param overall logical. If TRUE, then the EC index will be added to the result which is transformed into a list. Default equal to FALSE
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to ha).
#' @param rasterparallel logical. If nodes is "raster" then you can use this argument to assign the metrics values to the nodes raster. It is useful when raster resolution is less than 100 m2.
#' @param write character. Write output shapefile and overall table (if TRUE overall argument).
#'   It is necessary to specify the "Folder direction" + "Initial prefix",  for example, "C:/ejemplo".
#' @references Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. 2006. Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology 21 (7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.
#' Hanski, I. and Ovaskainen, O. 2000. The metapopulation capacity of a fragmented landscape. Nature 404: 755–758.
#' @export
#' @examples
#' library(Makurhini)
#' data("vegetation_patches", package = "Makurhini")
#' nrow(vegetation_patches) # Number of patches
#' #Two distance threshold,
#' IIC <- MK_dPCIIC(nodes = vegetation_patches, attribute = NULL,
#'                 distance = list(type = "centroid"),
#'                 metric = "IIC", distance_thresholds = c(5000, 10000)) # 5 and 10 km
#' IIC
#' #Using raster
#' data("raster_vegetation_patches", package = "Makurhini")
#' PC <- MK_dPCIIC(nodes = raster_vegetation_patches, attribute = NULL,
#'                 distance = list(type = "centroid"),
#'                 metric = "PC", probability = 0.5,
#'                 overall = TRUE,
#'                 distance_thresholds = 40000) # 40 km
#' PC
#'
#' @importFrom dplyr progress_estimated
#' @importFrom methods as
#' @importFrom utils write.table
#' @importFrom purrr map map_dbl
#' @importFrom raster as.matrix extent raster stack extent<- writeRaster reclassify crs crs<- unique
#' @importFrom future multiprocess plan
#' @importFrom furrr future_map
#' @importFrom igraph graph.adjacency shortest.paths E
#' @importFrom data.table data.table

MK_dPCIIC <- function(nodes, attribute  = NULL,
                      area_unit = "ha", restauration = NULL,
                      distance = list(type= "centroid", resistance = NULL),
                      metric = c("IIC", "PC"),
                      probability = NULL, distance_thresholds = NULL,
                      overall = FALSE, onlyoverall=FALSE,
                      LA = NULL, rasterparallel = FALSE, write = NULL) {
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (!metric %in% c("IIC", "PC")) {
    stop("Type must be either 'IIC', or 'PC'")
  }

  if (isTRUE(unique(metric == c("IIC", "PC")))) {
    metric = "IIC"
  }

  if (metric == "PC") {
    if (!is.null(probability) & !is.numeric(probability)) {
      stop("error missing probability")
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
  if(class(nodes)[1] == "SpatialPolygonsDataFrame" | class(nodes)[1] == "sf"){
    if (nrow(nodes) < 2) {
      stop("error, you need more than 2 nodes")
    }
    if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
      nodes <- st_as_sf(nodes) %>% st_zm()
    }
    nodes$IdTemp <- 1:nrow(nodes)
    idT <- "IdTemp"
  }

  if(class(nodes)[1] == "RasterLayer"){
    if(length(raster::unique(nodes)) < 2){
      stop("error, you need more than 2 nodes")
    }
    idT <- NULL
  }

  if(class(nodes)[1] == "data.frame"){
    attribute_1 <- nodes

    if(nrow(nodes) != nrow(attribute_1)){
      stop("nrow(nodes) != nrow(attribute)")
    }
    idT <- NULL

  } else {

    attribute_1 <- nodesfile(nodes, id = idT, attribute = attribute, area_unit = area_unit,
                             multiple = NULL, restauration = restauration,
                             prefix=NULL, write = NULL)
    }

    if(!is.null(restauration)){
    if(isFALSE(ncol(attribute_1) == 3 & isTRUE(unique(unique(attribute_1[,3])<=1)))){
      stop("error, review restauration argument, unique values must be 0 and 1")
    }
  }

  if(is.matrix(distance)){
    dist <- distance
    if(length(rownames(dist)) != nrow(nodes)){
      stop("nrow(dist) != nrow(nodes)")
    }

    if(isFALSE(unique(rownames(dist) %in% nodes[[1]]))){
      stop("Distance matrix and nodes have different id")
    }

  } else {
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
  }

  if(!is.null(idT)){
    rownames(dist) <- nodes$IdTemp
    colnames(dist) <- nodes$IdTemp
  }

  pb <- progress_estimated(length(distance_thresholds), 0)
  result <- map(distance_thresholds, function(x){

    if(length(distance_thresholds)>1){
      pb$tick()$print()
    }

    nodes.2 <- nodes

    if(metric == "IIC"){
      Adj_matr <- dist*0
      Adj_matr[dist < x] <- 1
      diag(Adj_matr) <- 0
      graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected"), error = function(err) err)

      if (inherits(graph_nodes, "error")) {
        stop("graph adjacency error")
      }

      nl.mat <- tryCatch(shortest.paths(graph_nodes), error = function(err) err)

      if (inherits(nl.mat, "error")) {
        stop("graph shortest.paths error")
      }

      nl.mat[is.infinite(nl.mat)] <- 1000

      attribute_2 <- attribute_1[,2]
      IICmat <- outer(attribute_2, attribute_2) / (1 + nl.mat)

      #p1
      IICnum <- sum(IICmat)

      #p2
      if(isTRUE(overall)){
        if(!is.null(LA)){
          if(LA > sum(attribute_2)){
            IIC <- IICnum / (LA^2)
            IIC_EC <- sqrt(IICnum)
            overall.2 = data.table(Index = c("IICnum", "EC(IIC)", "IIC"),
                                               Value = c(IICnum, IIC_EC, IIC))
          } else {
            stop("LA must be greater than the sum of the attributes of the nodes")
          }
        } else {
          IIC_EC <- sqrt(IICnum)
          overall.2 = data.table(Index = c("IICnum", "EC(IIC)"), Value = c(IICnum, IIC_EC))
        }
      }

      #p3
      if(isFALSE(onlyoverall)){
        N <- nrow(Adj_matr)
        dIIC <- map_dbl(1:N, function(i){
          mat.i <- Adj_matr[-i,-i]
          attribute.i <- attribute_2[-i]
          IIC.g.i <- graph.adjacency(mat.i, mode = "undirected")
          nl.mat.i <- shortest.paths(IIC.g.i)
          IICmat.i <- outer(attribute.i, attribute.i) / (1 + nl.mat.i)
          IICnum.i <- sum(IICmat.i)
          dIIC <- (IICnum-IICnum.i) / IICnum * 100
          return(dIIC)})

        #p4
        dIICintra <- attribute_2^2 / sum(IICmat) * 100
        dIICflux <- 2*(rowSums(IICmat) - attribute_2^2)/sum(IICmat) * 100
        dIICconnnector <- map_dbl(dIIC - dIICintra - dIICflux, function(x){if(x<0){0}else{x}})
        metric_conn <- as.data.frame(cbind(attribute_1[,1], dIIC, dIICintra, dIICflux, dIICconnnector))
        names(metric_conn)[1] <- c("Id")

        if(!is.null(restauration)){
          #T1 = 0
          IdT0 <- attribute_1[which(attribute_1$rest == 0),1]
          nl.mat.T0 <- nl.mat[which(rownames(nl.mat) %in% IdT0), which(colnames(nl.mat) %in% IdT0)]

          attribute_2 <- attribute_1[which(attribute_1[,1] %in% IdT0), 2]
          IICmat.T0 <- outer(attribute_2, attribute_2) / (1 + nl.mat.T0)
          #p1
          IICnum.T0 <- sum(IICmat.T0)

          #p3
          `%notin%` <- Negate(`%in%`)
          N <- which(attribute_1$rest == 1)
          attribute_2 <- attribute_1[, 2]
          dIICres <- rep(0, nrow(attribute_1))

          dIICres[N] <- map_dbl(N, function(i){
            mat.i <- Adj_matr[- N[which(N %notin% i)],-N[which(N %notin% i)]]
            attribute.i <- attribute_2[-N[which(N %notin% i)]]
            IIC.g.i <- graph.adjacency(mat.i, mode = "undirected")
            nl.mat.i <- shortest.paths(IIC.g.i)
            IICmat.i <- outer(attribute.i, attribute.i) / (1 + nl.mat.i)
            IICnum.i <- sum(IICmat.i)
            dIIC <- (IICnum.i- IICnum.T0) / IICnum.T0 * 100
            return(dIIC)})
          metric_conn <- cbind(metric_conn, dIICres)

        }

        if(!is.null(idT)){
          nodes.2$dIIC <- dIIC
          nodes.2$dIICintra <- dIICintra
          nodes.2$dIICflux <- dIICflux
          nodes.2$dIICconnnector <- dIICconnnector
          if(!is.null(restauration)){
            nodes.2$dIICres <- dIICres
          }
          nodes.2 <- nodes.2[moveme(names(nodes.2), "geometry last")]
          nodes.2$IdTemp <- NULL

          if(!is.null(write)){
            write_sf(nodes.2, paste0(write, "_", "d", x,  ".shp"), delete_layer = TRUE)
          }

        } else {
          if(class(nodes) == "data.frame"){

          } else{
          rp <- raster::unique(nodes)
          rp <- as.vector(rp)
          rp <- rp[which(!is.na(rp))]

          if(isTRUE(rasterparallel)){
            m <- matrix(nrow = nrow(attribute_1), ncol = 2)
            m[,1] <- attribute_1[,1]

            plan(strategy = multiprocess)
            r_metric <- future_map(2:ncol(metric_conn), function(c){
              x1 <- metric_conn[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)}, .progress = TRUE)

          } else {
            m <- matrix(nrow = nrow(attribute_1), ncol = 2)
            m[,1] <- attribute_1[,1]
            r_metric <- map(2:ncol(metric_conn), function(c){
              x1 <- metric_conn[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
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
        }
      }

    } else if (metric == "PC"){
      Adj_matr <- dist*0
      #negative kernel density
      if(is.null(probability)){
        k = (1 / x)
        Adj_matr <- exp(-k * dist)
      } else {
        k = log(probability)/x
        Adj_matr <- exp(k * dist)
      }
      diag(Adj_matr) <- 0

      #adjacency
      graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected", weighted = TRUE), error = function(err) err)

      if (inherits(graph_nodes, "error")) {
        stop("graph adjacency error")
      }

      #product of shortest paths
      pij.mat <- tryCatch(shortest.paths(graph_nodes, weights = -log(E(graph_nodes)$weight)), error = function(err) err)

      if (inherits(pij.mat, "error")) {
        stop("graph shortest.paths error")
      } else {
        pij.mat <- exp(-pij.mat)
      }

      attribute_2 <- attribute_1[,2]
      PCmat <- outer(attribute_2, attribute_2) * pij.mat

      #p1
      PCnum <- sum(PCmat)

      #p2
      if(isTRUE(overall)){
        if(!is.null(LA)){
          if(LA > sum(attribute_2)){
            PC <- PCnum / (LA^2)
            PC_EC <- sqrt(PCnum)
            overall.2 = data.table(Index = c("PCnum", "EC(PC)", "PC"),
                                               Value = c(PCnum, PC_EC, PC))
          } else {
            stop("LA must be greater than the sum of the attributes of the nodes")
          }
        } else {
          PC_EC <- sqrt(PCnum)
          overall.2 = data.table(Index = c("PCnum", "EC(PC)"), Value = c(PCnum, PC_EC))
        }
      }

      if(isFALSE(onlyoverall)){
        #p3
        N <- nrow(Adj_matr)
        dPC <- map_dbl(1:N, function(i){
          mat.i <- Adj_matr[-i,-i]
          attribute.i <- attribute_2[-i]
          PC.i <- graph.adjacency(mat.i, mode = "undirected", weighted = TRUE)
          pij.mat.i <- shortest.paths(PC.i, weights = -log(E(PC.i)$weight))
          pij.mat.i <- exp(-pij.mat.i)
          PCmat.i <- outer(attribute.i, attribute.i) * pij.mat.i
          PCnum.i <- sum(PCmat.i)
          dPC <- (PCnum-PCnum.i) / PCnum * 100
          return(dPC)})

        #p4
        dPCintra <- attribute_2^2 / sum(PCmat) * 100
        dPCflux <- 2*(rowSums(PCmat) - attribute_2^2)/sum(PCmat) * 100
        dPCconnector <- map_dbl(dPC - dPCintra - dPCflux, function(x){if(x<0){0}else{x}})
        metric_conn <- as.data.frame(cbind(attribute_1[,1], dPC, dPCintra, dPCflux, dPCconnector))
        names(metric_conn)[1] <- "Id"

        if(!is.null(restauration)){
          #T1 = 0
          IdT0 <- attribute_1[which(attribute_1$rest == 0), 1]
          pij.mat.T0 <- pij.mat[which(rownames(pij.mat) %in% IdT0), which(colnames(pij.mat) %in% IdT0)]

          attribute_2 <- attribute_1[which(attribute_1[,1] %in% IdT0), 2]
          PCmat.T0 <- outer(attribute_2, attribute_2) * pij.mat.T0
          #p1
          PCnum.T0 <- sum(PCmat.T0)

          #p3
          `%notin%` <- Negate(`%in%`)
          N <- which(attribute_1$rest == 1)
          attribute_2 <- attribute_1[, 2]
          dPCres <- rep(0, nrow(attribute_1))
          dPCres[N] <- map_dbl(N, function(i){
            mat.i <- Adj_matr[- N[which(N %notin% i)],-N[which(N %notin% i)]]
            attribute.i <- attribute_2[-N[which(N %notin% i)]]
            PC.i <- graph.adjacency(mat.i, mode = "undirected", weighted = TRUE)
            pij.mat.i <- shortest.paths(PC.i, weights = -log(E(PC.i)$weight))
            pij.mat.i <- exp(-pij.mat.i)
            PCmat.i <- outer(attribute.i, attribute.i) * pij.mat.i
            PCnum.i <- sum(PCmat.i)
            dPC <- (PCnum.i- PCnum.T0) / PCnum.T0 * 100
            return(dPC)})
          metric_conn <- cbind(metric_conn, dPCres)
        }

        if(!is.null(idT)){
          nodes.2$dPC <- dPC
          nodes.2$dPCintra <- dPCintra
          nodes.2$dPCflux <- dPCflux
          nodes.2$dPCconnector <- dPCconnector
          if(!is.null(restauration)){
            nodes.2$dPCres <- dPCres
          }
          nodes.2 <- nodes.2[moveme(names(nodes.2), "geometry last")]
          nodes.2$IdTemp <- NULL
          if(!is.null(write)){
            write_sf(nodes.2, paste0(write, "_", "d", x,  ".shp"), delete_layer = TRUE)
          }
        } else {
          if(class(nodes) == "data.frame"){
           nodes.2 <- cbind(nodes, metric_conn)
          } else {
            rp <- raster::unique(nodes)
            rp <- as.vector(rp)
            rp <- rp[which(!is.na(rp))]

            if(isTRUE(rasterparallel)){
              m <- matrix(nrow = nrow(attribute_1), ncol = 2)
              m[,1] <- attribute_1[,1]

              plan(strategy = multiprocess)
              r_metric <- future_map(2:ncol(metric_conn), function(c){
                x1 <- metric_conn[,c(1, c)]
                for(i in rp){
                  m[which(m == i),2] <- x1[which(x1[,1]== i),2]
                }
                x1 <- reclassify(nodes, rcl = m)
                return(x1)}, .progress = TRUE)

            } else {
              m <- matrix(nrow = nrow(attribute_1), ncol = 2)
              m[,1] <- attribute_1[,1]
              r_metric <- map(2:ncol(metric_conn), function(c){
                x1 <- metric_conn[,c(1, c)]
                for(i in rp){
                  m[which(m == i),2] <- x1[which(x1[,1]== i),2]
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
        }
      }

    } else {
      stop("error, select a correct metric")
    }

    if(isFALSE(onlyoverall)){
      if(isTRUE(overall)){
        result_metric <- list(nodes.2, overall.2)
        names(result_metric) <- c(paste0("node_importances_d",x), paste0("overall_d", x))
      } else {
        result_metric <- nodes.2
      }
    } else {
      result_metric <- overall.2
    }
    return(result_metric) })

  if(isFALSE(onlyoverall)){
  if (isTRUE(isFALSE(overall) && length(distance_thresholds) == 1)){
    result <- result[[1]]

  } else {
    if(length(distance_thresholds) == 1){
      result <- result[[1]]
    } else {
      names(result) <- paste0("d", distance_thresholds)
    }
  }
  } else {
    if(length(distance_thresholds) == 1){
      result <- result[[1]]
    } else {
      names(result) <- paste0("d", distance_thresholds)
    }

    }

  return(result)
}
