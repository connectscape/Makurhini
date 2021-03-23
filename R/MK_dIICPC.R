#' Estimate the integral index of connectivity (IIC) or the probability of connectivity (PC)
#'
#' Use this function to calculate the PC and IIC indexes under one or several distance thresholds.
#' @param nodes Nodes, patches, fragments, etc. Object of class sf, SpatialPolygonsDataFrame, raster or data frame.
#' If spatial layer, then it must be in a projected coordinate system. If nodes is a raster layer then raster values (Integer)
#' will be taken as "id". If nodes is a data frame, then it must have at least two columns, first column with nodes "id" and second with the attribute.
#' If you use the restoration argument, then you must have a dataframe with three columns where the third column is equal the restoration values
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#'   If nodes is a data frame then it must have two columns where second column is the attribute.
#' @param area_unit character. If attribute is NULL you can set an area unit (e.g., "km2", "cm2", "ha";
#'  see Makurhini::unit_convert). Default equal to hectares "m2".
#' @param restoration character or vector. If nodes is a shappefile then you must specify the name of the column
#' with restoration value. If nodes is a raster layer then must be a numeric vector with restoration values
#' to each node in the raster. Binary values (0,1), where 1 = existing nodes in the landscape, and 0 = a new node
#'  to add to the initial landscape (restored).
#'  If nodes is a data frame then it must have three columns where third column has restoration values.
#' @param distance list or distance matrix. If it is a list, then it must be a set of distance parameters. For example: type, resistance,or keep. For "type" choose one
#' of the distances: "centroid" (faster), "edge", "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have
#' to use the "resistance" argument. To see more arguments see the ?distancefile.
#' If it is a distance matrix, then the number of columns and rows must be equal the number of nodes and "id".
#' @param metric character. Choose a connectivity metric: "IIC" considering topological distances or "PC" considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' If probability = NULL, then it will be the inverse of the mean dispersal distance for the species (1/α; Hanski and Ovaskainen 2000).
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000). If NULL then the mean
#'  distance between nodes will be estimated and used.
#' @param overall logical. If TRUE, then the EC index will be added to the result which is transformed into a list. Default equal to FALSE
#' @param onlyoverall logical. If TRUE, then only overall metrics will be calculated.
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to m2).
#' @param rasterparallel numeric. Specify the number of cores to use for parallel processing, default = NULL.
#' If nodes is "raster" then you can use this argument to assign the metrics values to the nodes
#'  raster. It is useful when raster resolution is less than 100 m2.
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
#'                 area_unit = "m2",
#'                 distance = list(type = "centroid"),
#'                 metric = "IIC", distance_thresholds = c(5000, 10000)) # 5 and 10 km
#' IIC
#' plot(IIC$d5000["dIIC"], breaks = "jenks")
#' plot(IIC$d5000["dIICflux"], breaks = "jenks")
#' plot(IIC$d5000["dIICconnector"])
#'
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
#' @importFrom future multiprocess plan availableCores
#' @importFrom furrr future_map
#' @importFrom igraph graph.adjacency shortest.paths E
#' @importFrom data.table data.table fwrite
#' @importFrom sf write_sf st_as_sf st_zm

MK_dPCIIC <- function(nodes, attribute  = NULL,
                      area_unit = "m2", restoration = NULL,
                      distance = list(type= "centroid", resistance = NULL),
                      metric = c("IIC", "PC"),
                      probability = NULL, distance_thresholds = NULL,
                      overall = FALSE, onlyoverall=FALSE,
                      LA = NULL, rasterparallel = NULL, write = NULL) {
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
      nodes <- st_as_sf(nodes)
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

    if(dim(nodes)[2] < 2){
      stop("You need two or three columns, see the ?MK_dPCIIC")
    }

    idT <- NULL

  } else {
    attribute_1 <- nodesfile(nodes, id = idT, attribute = attribute,
                             area_unit = area_unit,
                             restoration = restoration,
                             write = NULL)
  }

  if(!is.null(restoration)){
    if(isFALSE(ncol(attribute_1) == 3 &
               isTRUE(unique(unique(attribute_1[, 3]) <= 1)))){
      stop("error, review restoration argument, unique values must be 0 and 1")
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
    dist <- distancefile(nodes = nodes,  id = idT,
                         type = distance$type,
                         distance_unit = distance$distance_unit,
                         keep = distance$keep,
                         resistance = distance$resistance,
                         resist.units = distance$resist.units,
                         CostFun = distance$CostFun,
                         ngh = distance$ngh,
                         mask = distance$mask,
                         threshold = distance$threshold,
                         geometry_out = distance$geometry_out,
                         bounding_circles = distance$bounding_circles,
                         parallel = distance$parallel,
                         edgeParallel = distance$edgeParallel,
                         least_cost.java = distance$least_cost.java,
                         cores.java = distance$cores.java, ram.java = distance$ram.java,
                         pairwise = FALSE,
                         write = NULL)
  }

  pb <- progress_estimated(length(distance_thresholds), 0)

  result <- map(distance_thresholds, function(x){

    if(length(distance_thresholds)>1){
      pb$tick()$print()
    }

    nodes.2 <- nodes
    Adj_matr <- dist * 0

    if(is.null(x)){
      x <- mean(dist)
    }

    if(metric == "IIC"){
      Adj_matr[dist < x] <- 1

    } else {
      #negative kernel density
      if(is.null(probability)){
        k = (1 / x)
        Adj_matr <- exp(-k * dist)
      } else {
        Adj_matr <- exp((dist * log(probability))/x)
      }
    }

    diag(Adj_matr) <- 0

    mode(Adj_matr) <- "numeric"

    #adjacency
    graph_nodes <- tryCatch(graph.adjacency(Adj_matr, mode = "undirected",
                                            weighted = if(metric == "IIC"){NULL}else{TRUE}),
                            error = function(err) err)

    if (inherits(graph_nodes, "error")) {
      stop("graph adjacency error")
    }

    mat1 <- tryCatch(shortest.paths(graph_nodes,
                                    weights = if(metric == "IIC"){NULL
                                    } else {-log(E(graph_nodes)$weight)}),
                     error = function(err) err)

    if (inherits(mat1, "error")) {
      stop("graph shortest.paths error")
    } else {
      if(metric == "PC"){
        mat1 <- exp(-mat1)
      }
    }

    mat1[is.infinite(mat1)] <- 1000

    attribute_2 <- attribute_1[,2]

    if(metric == "IIC"){
      mat2 <- outer(attribute_2, attribute_2) / (1 + mat1)
    } else {
      mat2 <- outer(attribute_2, attribute_2) * mat1
    }

    #p1
    num <- sum(mat2)

    #p2
    last <- if(metric == "IIC"){"IIC"} else {"PC"}
    if(isTRUE(overall)){
      if(!is.null(LA)){
        if(LA > sum(attribute_2)){
          Ind <- num / (LA^2)
          EC <- sqrt(num)
          overall.2 = data.table(Index = c(paste0(last, "num"),
                                           paste0("EC(", last, ")"),
                                           last),
                                 Value = c(num, EC, Ind))
        } else {
          stop("LA must be greater than the sum of the attributes of the nodes")
        }
      } else {
        EC <- sqrt(num)
        overall.2 = data.table(Index = c(paste0(last, "num"),
                                         paste0("EC(", last, ")")),
                               Value = c(num, EC))
      }
    }

    #
    if(isFALSE(onlyoverall)){
      #p3
      N <- nrow(Adj_matr)
      delta <- map_dbl(1:N, function(i){
        mat.i <- Adj_matr[-i,-i]
        attribute.i <- attribute_2[-i]
        g.i <- graph.adjacency(mat.i, mode = "undirected",
                               weighted = if(metric == "IIC"){NULL}else{TRUE})
        mat.i <- shortest.paths(g.i, weights = if(metric == "IIC"){NULL}else{-log(E(g.i)$weight)})
        if(metric == "PC"){
          mat.i <- exp(-mat.i)
        }
        if(metric == "IIC"){
          mat.i <- outer(attribute.i, attribute.i) / (1 + mat.i)
        } else {
          mat.i <- outer(attribute.i, attribute.i) * mat.i
        }
        num.i <- sum(mat.i)
        dC <- (num - num.i) / num * 100
        return(dC)})

      #p4
      dintra <- attribute_2^2 / num * 100
      dflux <- 2 * (rowSums(mat2) - attribute_2^2)/ num * 100
      dconnector <- map_dbl(delta - dintra - dflux, function(x){if(x < 0){0} else {x}})
      metric_conn <- as.data.frame(cbind(attribute_1[,1], delta, dintra, dflux, dconnector))
      names(metric_conn)[1] <- c("Id")

      #p5
      if(!is.null(restoration)){
        #T1 = 0
        IdT0 <- attribute_1[which(attribute_1[,3] == 0), 1]
        mat.T0 <- mat1[which(rownames(mat1) %in% IdT0), which(colnames(mat1) %in% IdT0)]

        attribute_2 <- attribute_1[which(attribute_1[,1] %in% IdT0), 2]

        if(metric == "IIC"){
          mat.T0 <- outer(attribute_2, attribute_2) / (1 + mat.T0)
        } else {
          mat.T0 <- outer(attribute_2, attribute_2) * mat.T0
        }

        #p1
        num.T0 <- sum(mat.T0)

        #p3
        `%notin%` <- Negate(`%in%`)
        N <- which(attribute_1[,3] == 1)
        attribute_2 <- attribute_1[, 2]
        dres <- rep(0, nrow(attribute_1))

        dres[N] <- map_dbl(N, function(i){
          mat.i <- Adj_matr[- N[which(N %notin% i)],-N[which(N %notin% i)]]
          attribute.i <- attribute_2[-N[which(N %notin% i)]]
          g.i <- graph.adjacency(mat.i, mode = "undirected",
                                 weighted = if(metric == "IIC"){NULL}else{TRUE})
          mat.i <- shortest.paths(g.i, weights = if(metric == "IIC"){NULL}else{-log(E(g.i)$weight)})
          if(metric == "PC"){
            mat.i <- exp(-mat.i)
          }
          if(metric == "IIC"){
            mat.i <- outer(attribute.i, attribute.i) / (1 + mat.i)
          } else {
            mat.i <- outer(attribute.i, attribute.i) * mat.i
          }

          num.i <- sum(mat.i)
          dC <- (num.i- num.T0) / num.T0 * 100
          return(dC)})
        metric_conn <- cbind(metric_conn, dres)
      }

      if(!is.null(idT)){
        nodes.2 <- cbind(nodes.2, metric_conn)

        names(nodes.2)[which(names(nodes.2) %in%
                               c("delta", "dintra", "dflux", "dconnector"))] <- c(paste0("d", metric),
                                                                                  paste0("d", metric, "intra"),
                                                                                  paste0("d", metric, "flux"),
                                                                                  paste0("d", metric, "connector"))

        if(!is.null(restoration)){
          nodes.2$dres <- dres
          names(nodes.2)[which(names(nodes.2) == "dres")] <- paste0("d", metric, "res")
        }

        nodes.2$IdTemp <- NULL
        nodes.2$Id <- NULL
        nodes.2 <- nodes.2[,which(names(nodes.2) != "geometry")]

        if(!is.null(write)){
          write_sf(nodes.2, paste0(write, "_", "d", x,  ".shp"), delete_layer = TRUE)
        }

      } else {
        names(metric_conn)[2:ncol(metric_conn)] <- c(paste0("d", metric),
                                                     paste0("d", metric, "intra"),
                                                     paste0("d", metric, "flux"),
                                                     paste0("d", metric, "connector"))
        if(class(nodes) == "data.frame"){
          nodes.2 <- cbind(nodes, metric_conn)
          nodes.2$Id <- NULL
          nodes.2$IdTemp <- NULL

          if(!is.null(write)){
            fwrite(nodes.2, paste0(write, "_", "d", x,  ".csv"))
          }
        } else {
          rp <- raster::unique(nodes)
          rp <- as.vector(rp)
          rp <- rp[which(!is.na(rp))]

          if(!is.null(rasterparallel)){
            m <- matrix(nrow = nrow(attribute_1), ncol = 2)
            m[,1] <- attribute_1[,1]

            works <- as.numeric(availableCores())-1
            works <-  if(rasterparallel > works){works}else{rasterparallel}
            plan(strategy = multiprocess, gc = TRUE, workers = works)
            r_metric <- future_map(2:ncol(metric_conn), function(c){
              x1 <- metric_conn[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)}, .progress = TRUE)
            close_multiprocess(works)

          } else {
            m <- matrix(nrow = nrow(attribute_1), ncol = 2)
            m[,1] <- attribute_1[,1]
            r_metric <- lapply(2:ncol(metric_conn), function(c){
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

      if(isTRUE(overall)){
        result_metric <- list(nodes.2, overall.2)
        names(result_metric) <- c(paste0("node_importances_d",x), paste0("overall_d", x))
      } else {
        result_metric <- nodes.2
      }

    } else {
      result_metric <- overall.2
    }

    #
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
