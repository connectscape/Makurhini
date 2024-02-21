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
#'  see Makurhini::unit_convert). Default equal to square meters "m2".
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
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000). If NULL then the median
#'  distance between nodes will be estimated and used.Also, you can use the dispersal_distance() function to estimate the a dispersal distance using the species home range.
#' @param overall logical. If TRUE, then the EC index will be added to the result which is transformed into a list. Default equal to FALSE
#' @param onlyoverall logical. If TRUE, then only overall metrics will be calculated.
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to m2).
#' @param parallel numeric. Specifies the number of cores to parallelize the estimation of shortest paths.It is useful when you have more than 1,000 nodes.
#' @param rasterparallel numeric. Specify the number of cores to use for parallel processing, default = NULL.
#' If nodes is "raster" then you can use this argument to assign the metrics values to the nodes
#'  raster. It is useful when raster resolution is less than 100 m2.
#' @param write character. Write output shapefile and overall table (if TRUE overall argument).
#'   It is necessary to specify the "Folder direction" + "Initial prefix",  for example, "C:/ejemplo".
#' @param intern logical. Show the progress of the process, default = TRUE. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @references Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. 2006. Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology 21 (7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.
#' Hanski, I. and Ovaskainen, O. 2000. The metapopulation capacity of a fragmented landscape. Nature 404: 755–758.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' data("vegetation_patches", package = "Makurhini")
#' nrow(vegetation_patches) # Number of patches
#' #Two distance threshold,
#' IIC <- MK_dPCIIC(nodes = vegetation_patches, attribute = NULL,
#'                 area_unit = "m2",
#'                 distance = list(type = "centroid"),
#'                 parallel = NULL,
#'                 metric = "IIC", distance_thresholds = c(10000, 20000)) #10,20 km
#' IIC
#' plot(IIC$d20000["dIIC"], breaks = "jenks")
#' plot(IIC$d20000["dIICflux"], breaks = "jenks")
#' plot(IIC$d20000["dIICconnector"], breaks = "jenks")
#' #Using raster
#' data("raster_vegetation_patches", package = "Makurhini")
#' PC <- MK_dPCIIC(nodes = raster_vegetation_patches, attribute = NULL,
#'                 distance = list(type = "centroid"),
#'                 metric = "PC", probability = 0.5,
#'                 overall = TRUE,
#'                 distance_thresholds = 40000) # 40 km
#' PC$overall_d40000
#' plot(PC$node_importances_d40000)
#' }
#' @note Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom purrr map_dbl
#' @importFrom raster as.matrix extent raster stack extent<- writeRaster reclassify crs crs<- unique
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map future_map_dbl
#' @importFrom igraph graph.adjacency shortest.paths E as_ids V
#' @importFrom sf write_sf st_as_sf st_zm
#' @importFrom magrittr %>%
#' @importFrom stats median
MK_dPCIIC <- function(nodes, attribute  = NULL,
                      area_unit = "m2", restoration = NULL,
                      distance = list(type= "centroid", resistance = NULL),
                      metric = c("IIC", "PC"),
                      probability = NULL, distance_thresholds = NULL,
                      overall = FALSE, onlyoverall = FALSE,
                      LA = NULL,
                      parallel = NULL,
                      rasterparallel = NULL, write = NULL,
                      intern = TRUE) {
  . = NULL
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (isTRUE(unique(c("IIC", "PC") %in% metric))) {
    metric = "IIC"; instr_dist <- TRUE
  } else {
    if (metric == "PC") {
      if (!is.null(probability) & !is.numeric(probability)) {
        stop("error missing probability")
      }
      instr_dist <- FALSE
    } else if (metric == "IIC"){
      instr_dist <- TRUE
    } else {
      stop("Type must be either 'IIC', or 'PC'")
    }
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(isFALSE(parallel)){
    parallel <- NULL
  } else if (isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
  } else if((!is.null(parallel))){
    if(!is.numeric(parallel)){
      stop("if you use parallel argument then you need a numeric value")
    }
  } else {
    parallel <- NULL
  }

  #Nodes
  if(class(nodes)[1] == "SpatialPolygonsDataFrame" | class(nodes)[1] == "sf"){
    if (nrow(nodes) < 2) {
      stop("error, you need more than 2 nodes")
    }

    if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
      nodes <- st_as_sf(nodes)
    }

    nodes$IdTemp <- 1:nrow(nodes); idT <- "IdTemp"
    attribute_1 <- nodesfile(nodes, id = idT, attribute = attribute,
                             area_unit = area_unit,
                             restoration = restoration,
                             write = NULL)
  } else if (class(nodes)[1] == "RasterLayer"){
    if(length(raster::unique(nodes)) < 2){
      stop("error, you need more than 2 nodes")
    }
    idT <- NULL
    attribute_1 <- nodesfile(nodes, id = idT, attribute = attribute,
                             area_unit = area_unit,
                             restoration = restoration,
                             write = NULL)
  } else if (class(nodes)[1] == "data.frame"){
    attribute_1 <- nodes; idT <- NULL

    if(dim(nodes)[2] < 2){
      stop("You need two or three columns, see the ?MK_dPCIIC")
    }
  } else {
    stop("error missing file of nodes")
  }

  if(!is.null(restoration)){
    if(isFALSE(ncol(attribute_1) == 3 & isTRUE(unique(unique(attribute_1[, 3]) <= 1)))){
      stop("error, review restoration argument, unique values must be 0 and 1")
    }
  }

  #Distance
  if(is.matrix(distance)){
    dist <- distance

    if(unique(dim(dist)) != nrow(nodes)){
      stop("nrow(dist) != nrow(nodes)")
    }

    if(isFALSE(unique(rownames(dist) %in% nodes[[1]]))){
      stop("Distance matrix and nodes have different id")
    }

  } else {
    dist <- tryCatch(distancefile(nodes = nodes,
                         id = idT,
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
                         ActiveParallel = FALSE,
                         least_cost.java = distance$least_cost.java,
                         cores.java = distance$cores.java,
                         ram.java = distance$ram.java,
                         pairwise = FALSE,
                         write = NULL), error = function(err)err)
    #dist <- distance
    if(inherits(dist, "error")){
      close_multiprocess()
      stop(dist)
    }
  }

  if(is.null(distance_thresholds)){
    distance_thresholds <- median(dist)
  }

  if(length(distance_thresholds) > 1 & isTRUE(intern)){
    pb <- txtProgressBar(0,length(distance_thresholds), style = 3)
  }
  #x=1
  #result <- list(result_metric)
  result <- lapply(1:length(distance_thresholds), function(x){
    x.1 <- distance_thresholds[x]; nodes.2 <- nodes

    mat1 <- tryCatch(get_sdist(dist_nodes = dist,
                               attr_nodes = attribute_1[,2],
                               metric = metric,
                               probability = probability,
                               distance_threshold = x.1,
                               igraph_Dijkstra = instr_dist,
                               parallel = parallel,
                               loop = TRUE, G1 = 1000, G2 = 1000,
                               intern = intern), error = function(err)err)
    #mat1 <- smat

    if(inherits(mat1, "error")){
      stop("Error in short distance estimation")
    } else {
      attribute_2 <- attribute_1[,2]
      num <- sum(mat1)
    }

    # sum(2*(rowSums(mat2) - attribute_2^2)/ num * 100)
    # sum(((rowSums(mat2)+colSums(mat2)) - attribute_2^2))
    # last <- if(metric == "IIC"){"IIC"} else {"PC"}

    if(isTRUE(overall) | isTRUE(onlyoverall)){
      if(!is.null(LA)){
        if(LA > sum(attribute_2)){
          overall.2 <- data.frame(Index = c(paste0(metric, "num"),
                                            paste0("EC(", metric, ")"),
                                            metric),
                                  Value = c(num, sqrt(num), num / (LA^2)))
          if(overall.2$Value[3] > 1){
            overall.2$Value[3] <- 1
          }

          if(overall.2$Value[2] > LA){
            overall.2$Value[2] <- LA
          }

        } else {
          stop("LA must be greater than the sum of the attributes of the nodes")
        }
      } else {
        overall.2 <- data.frame(Index = c(paste0(metric, "num"),
                                          paste0("EC(", metric, ")")),
                                Value = c(num, sqrt(num)))
      }

      if (!is.null(write)){
        write.csv(overall.2, file = paste0(write, "_Overall","_", "d", x.1,  ".csv"))
      }
    }

    if(isFALSE(onlyoverall)){
      if(nrow(attribute_1) < 1000 & is.null(parallel)){
        delta <- map_dbl(1:nrow(attribute_1), function(i){
          attribute.i <- attribute_2[-i]
          mat.i <- tryCatch(get_sdist(dist_nodes = dist[-i,-i],
                                      attr_nodes = attribute.i,
                                      metric = metric,
                                      probability = probability,
                                      distance_threshold = x.1,
                                      igraph_Dijkstra = instr_dist,
                                      parallel = parallel, loop = TRUE,
                                      G1 = 1000, G2 = 1000, intern = FALSE), error = function(err)err)

          if(inherits(mat.i, "error")){
            stop("Error in short distance estimation")
          } else {
            num.i <- sum(mat.i); dC <- (num - num.i) / num * 100; return(dC)
          }
          })
      } else {
        #Probar parallel fuera de la funcion o dentro de
        if(!is.null(parallel)){
          works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}
          if(.Platform$OS.type == "unix") {
            strat <- future::multicore
          } else {
            strat <- future::multisession
          }
          plan(strategy = strat, gc = TRUE, workers = works)
          delta <- future_map_dbl(1:nrow(attribute_1), function(i){
            attribute.i <- attribute_2[-i]
            mat.i <- tryCatch(get_sdist(dist_nodes = dist[-i,-i],
                                        attr_nodes = attribute.i,
                                        metric = metric,
                                        probability = probability,
                                        distance_threshold = x.1,
                                        igraph_Dijkstra = instr_dist,
                                        parallel = NULL, loop = TRUE,
                                        G1 = 1000, G2 = 1000, intern = FALSE), error = function(err)err)

            if(inherits(mat.i, "error")){
              stop("Error in short distance estimation")
            } else {
              num.i <- sum(mat.i); dC <- (num - num.i) / num * 100; return(dC)
            }
            }, .progress = intern)
        } else {
          delta <- map_dbl(1:nrow(attribute_1), function(i){
            attribute.i <- attribute_2[-i]
            mat.i <- tryCatch(get_sdist(dist_nodes = dist[-i,-i],
                                        attr_nodes = attribute.i,
                                        metric = metric,
                                        probability = probability,
                                        distance_threshold = x.1,
                                        igraph_Dijkstra = instr_dist,
                                        parallel = parallel, loop = TRUE,
                                        G1 = 1000, G2 = 1000, intern = FALSE), error = function(err)err)

            if(inherits(mat1, "error")){
              stop("Error in short distance estimation")
            } else {
              num.i <- sum(mat.i); dC <- (num - num.i) / num * 100; return(dC)
            }
            })
        }
      }
      dintra <- round((((attribute_2^2)) / num) * 100, 7); dflux <- round(2*(rowSums(mat1) - attribute_2^2)/ num * 100, 7)
      dconnector <- round(map_dbl(delta - ((((attribute_2^2)) / num) * 100) - (2*(rowSums(mat1) - attribute_2^2)/ num * 100), function(y){if(y < 0){0} else {y}}), 10)
      metric_conn <- data.frame("IdTemp2" = attribute_1[,1], "delta" = round(delta, 7),
                                "dintra" = dintra, "dflux" = dflux, "dconnector" = dconnector,
                                check.names = FALSE)

      if(!is.null(restoration)){
        IdT0 <- attribute_1[which(attribute_1[,3] == 0), 1]
        mat.T0 <- mat1[which(rownames(mat1) %in% IdT0), which(colnames(mat1) %in% IdT0)]

        attribute_2 <- attribute_1[which(attribute_1[,1] %in% IdT0), 2]

        if(metric == "IIC"){
          mat2 <- outer(attribute_2, attribute_2) / (1 + mat.T0); mat2[which(mat.T0 == 1000000)] <- 0
        } else {
          mat2 <- outer(attribute_2, attribute_2) * mat.T0
        }

        num.T0 <- sum(mat2); `%notin%` <- Negate(`%in%`)
        N <- which(attribute_1[,3] == 1); attribute_2 <- attribute_1[, 2]; dres <- rep(0, nrow(attribute_1))

        if(!is.null(parallel)){
          dres[N] <- future_map_dbl(N, function(i){
            dist.i <- dist[-N[which(N %notin% i)],-N[which(N %notin% i)]]; attribute.i <- attribute_2[-N[which(N %notin% i)]]

            mat.i <- tryCatch(get_sdist(dist_nodes = dist.i,
                                        attr_nodes = attribute.i,
                                        metric = metric,
                                        probability = probability,
                                        distance_threshold = x.1,
                                        igraph_Dijkstra = FALSE,
                                        parallel = NULL, loop = TRUE,
                                        G1 = 1000, G2 = 1000, intern = FALSE), error = function(err)err)

            if(inherits(mat.i, "error")){
              stop("Error in short distance estimation")
            } else {
              num.i <- sum(mat.i); dC <- (num.i- num.T0) / num.T0 * 100; return(dC)
            }
            })
          close_multiprocess(works)
        } else {
          dres[N] <- map_dbl(N, function(i){
            dist.i <- dist[-N[which(N %notin% i)],-N[which(N %notin% i)]]; attribute.i <- attribute_2[-N[which(N %notin% i)]]

            mat.i <- tryCatch(get_sdist(dist_nodes = dist.i,
                                        attr_nodes = attribute.i,
                                        metric = metric,
                                        probability = probability,
                                        distance_threshold = x.1,
                                        igraph_Dijkstra = FALSE,
                                        parallel = parallel, loop = TRUE,
                                        G1 = 1000, G2 = 1000, intern = FALSE), error = function(err)err)

            if(inherits(mat1, "error")){
              stop("Error in short distance estimation")
            } else {
              num.i <- sum(mat.i); dC <- (num.i- num.T0) / num.T0 * 100; return(dC)
            }
            })
        }

        metric_conn <- cbind(metric_conn, dres)
      } else {
        if(!is.null(parallel)){
          close_multiprocess(works)
        }
      }

      if(!is.null(idT)){
        nodes.2 <- cbind(nodes.2, metric_conn)
        names(nodes.2)[which(names(nodes.2) %in%
                               c("delta", "dintra", "dflux", "dconnector"))] <- c(paste0("d", metric),
                                                                                  paste0("d", metric, "intra"),
                                                                                  paste0("d", metric, "flux"),
                                                                                  paste0("d", metric, "connector"))

        if(!is.null(restoration)){
          names(nodes.2)[which(names(nodes.2) == "dres")] <- paste0("d", metric, "res")
        }

        nodes.2$IdTemp <- NULL; nodes.2$IdTemp2 <- NULL; nodes.2 <- nodes.2[,which(names(nodes.2) != "geometry")]

        if(!is.null(write)){
          write_sf(nodes.2, paste0(write, "_", "d", x.1,  ".shp"), delete_layer = TRUE)
        }

      } else {
        names(metric_conn)[2:ncol(metric_conn)] <- c(paste0("d", metric),
                                                     paste0("d", metric, "intra"),
                                                     paste0("d", metric, "flux"),
                                                     paste0("d", metric, "connector"))
        if(class(nodes)[1] == "data.frame"){
          nodes.2 <- cbind(nodes, metric_conn); nodes.2$IdTemp <- NULL; nodes.2$IdTemp2 <- NULL

          if(!is.null(write)){
            write.csv(nodes.2, paste0(write, "_", "d", x.1,  ".csv"))
          }
        } else {
          rp <- raster::unique(nodes); rp <- as.vector(rp); rp <- rp[which(!is.na(rp))]

          if(!is.null(rasterparallel)){
            m <- matrix(nrow = nrow(attribute_1), ncol = 2); m[,1] <- attribute_1[,1]

            works <- as.numeric(availableCores())-1; works <-  if(rasterparallel > works){works}else{rasterparallel}
            if(.Platform$OS.type == "unix") {
              strat <- future::multicore
            } else {
              strat <- future::multisession
            }
            plan(strategy = strat, gc = TRUE, workers = works)

            r_metric <- tryCatch(future_map(2:ncol(metric_conn), function(c){
              x1 <- metric_conn[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)}, .progress = TRUE), error = function(err) err)
            close_multiprocess(works)

          } else {
            m <- matrix(nrow = nrow(attribute_1), ncol = 2); m[,1] <- attribute_1[,1]
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

          nodes.2[[1]] <- nodes; nodes.2 <- stack(nodes.2)
          names(nodes.2) <- names(metric_conn); names(nodes.2)[1] <- "Id"

          if (!is.null(write)){
            n <- names(metric_conn); n <- lapply(as.list(2:length(n)), function(w){
              x1 <- nodes.2[[w]]; crs(x1) <- crs(nodes.2)
              writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x.1, ".tif"),
                          overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
            })
          }
        }
      }

      if(isTRUE(overall)){
        result_metric <- list(nodes.2, overall.2)
        names(result_metric) <- c(paste0("node_importances_d",x.1), paste0("overall_d", x.1))
      } else {
        result_metric <- nodes.2
      }

    } else {
      result_metric <- overall.2
    }

    if(length(distance_thresholds) > 1 & isTRUE(intern)){
      setTxtProgressBar(pb, x)
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
