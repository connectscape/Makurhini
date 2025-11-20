#' @title Estimate the probability of connectivity (PC) or the integral index of connectivity (IIC) and prioritize nodes
#'
#' @description This function calculates both the overall landscape connectivity, and the importance (contribution) of each node (or habitat patch) for maintaining landscape connectivity. It uses the PC and IIC indexes under one or several distance thresholds.
#' @details This function calculates the \bold{overall connectivity} and the \bold{importance or contribution of each node} to the overall landscape connectivity. The overall connectivity is computed using either the PC (Probability of Connectivity) or the IIC (Integral Index of Connectivity) metrics.
#' @param nodes Object containing nodes (e.g., habitat patches or fragments) information. It can be of the following classes:\cr
#' -   \code{Data.frame} with at least two columns: the first for node IDs and the second for attributes. If the `restoration` argument is used, the data frame must include a third column for restoration values.\cr
#' -   Spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.\cr
#' -   Raster (class \code{RasterLayer, SpatRaster}). It must be in a projected coordinate system. The values must be integers representing the ID of each habitat patch or node, with non-habitat areas represented by NA values (see \link[raster]{clump} or \link[terra]{patches}).
#' @param attribute \code{Character} or \code{vector}. If \code{NULL} (only applicable when \code{nodes} is of spatial data of vector or raster type) the area of the nodes will be used as the node attribute. The unit of area can be selected using the \code{area_unit} parameter. To use an alternative attribute, consider the class type of the object in the \code{nodes} parameter: \cr
#' -   If \code{nodes} is a spatial vector or data.frame, specify \bold{the name of the column} containing the attribute for the nodes. \cr
#' -   If \code{nodes} is a raster layer then it must be a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes. If the parameter \bold{weighted} is \code{TRUE} then the numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#' @param weighted \code{logical}. If the \code{nodes} are of raster type, you can weight the estimated area of each node by the attribute. When using this parameter the \code{attribute} parameter, which must be a vector of length equal to the number of nodes, usually has values between 0 and 1.
#' @param LA \code{numeric}. (\emph{optional, default = } \code{NULL}). The maximum landscape attribute, which is the attribute value that would correspond to a hypothetical habitat patch covering all the landscape with the best possible habitat, in which IIC and PC would be equal to 1. For example, if nodes attribute corresponds to the node area, then LA equals total landscape area. If nodes attribute correspond to a quality-weighted area and the quality factor ranges from 0 to 100, LA will be equal to 100 multiplied by total landscape area. The value of LA does not affect at all the importance of the nodes and is only used to calculate the overall landscape connectivity. If no LA value is entered (default) and  \code{overall = TRUE} or \code{onlyoverall = TRUE}, the function will only calculate the numerator of the global connectivity indices and the equivalent connected ECA or EC index.
#' @param area_unit \code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.
#' @param restoration \code{Character or vector}, (\emph{optional}). This parameter specifies the binary restoration value indicating whether each node is existing or hypothetically added for restoration: 1 for existing nodes in the landscape and 0 for nodes to be added to the initial landscape. If \code{NULL} (default), all nodes are considered existing (as if all restoration values were 1). Otherwise: \cr
#' -   If \code{nodes} is a shapefile (spatial vector) or a data.frame, specify the name of the column containing the restoration values. If `nodes` is a data.frame, it must have three columns: the first for node IDs, the second for attributes, and the third for restoration values.\cr
#' -   If \code{nodes} is a raster layer, provide a numeric vector with restoration values for each node in the raster.
#' @param onlyrestor \code{logical}. If \code{TRUE}, then only restoration metric will be calculated.
#' @param distance  A \code{matrix} or \code{list} to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. \cr
#' - If it is a \code{matrix}, then the number of columns and rows must be equal to the number of nodes. This distance matrix could be generated by the \link[Makurhini]{distancefile} function.\cr
#' - If it is a \code{list} of parameters, then it must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. For example: \code{distance(type = "least-cost", resistance = raster_resistance)}. \cr
#' To see more arguments see the \link[Makurhini]{distancefile} function.
#' @param metric A \code{character} indicating the connectivity metric to use: \code{"PC"} (the default and recommended) to calculate the probability of connectivity index, and \code{"IIC"} to calculate the binary integral index of connectivity.
#' @param probability A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the \code{"PC"} metric. If \code{probability = NULL}, then a probability of 0.5 will be used.
#' @param distance_thresholds A \code{numeric} indicating the dispersal distance or distances (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range.
#' @param threshold \code{numeric}. Pairs of nodes with a distance value greater than this threshold will be discarded in the analysis which can speed up processing. Can be the same length as the \code{distance_thresholds} parameter.
#' @param overall \code{logical}. If \code{TRUE}, then the overall metrics will be added to the result which is transformed into a list. Default equal to FALSE
#' @param onlyoverall \code{logical}. If \code{TRUE}, then only overall metrics will be calculated.
#' @param parallel  (\emph{optional, default =} \code{NULL}).
#' A \code{numeric} specifying the number of cores to parallelize the index estimation of the PC or IIC index and its deltas.Particularly useful when you have more than 1000 nodes. By default the analyses are not parallelized.
#' @param parallel_mode (\emph{optional, default =} \code{1}).
#' A \code{numeric} indicating the mode of parallelization: Mode \bold{1} (the default option, and recommended for less than 1000 nodes) parallelizes on the connectivity delta estimate, while Mode \bold{2} (recommended for more than 1000 nodes)
#' parallelizes on the shortest paths between vertices estimate.
#' @param write \code{Character} indicating the path and initial prefix of the objects to save, for example, \code{"C:/example/test_PC_"}. By default, nothing is saved. The saved objects are: \cr
#' -   The importance of each node or habitat patch. The format of the output file depends on the class of \code{nodes} (\code{shapefile, raster, or table}).\cr
#' -   Overall landscape connectivity table (if the \code{overall} argument is \code{TRUE}).
#' @param intern \code{logical}. Show the progress of the process, \code{default = TRUE}. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @param id_sel Internal use only, not for users.
#' @references
#' -   Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#' -   Pascual-Hortal, L. & Saura, S. 2006. Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology 21 (7): 959-967.\cr
#' -   Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.\cr
#' -   Saura, S., Bodin, Ö., & Fortin, M.-J. (2014). EDITOR’S CHOICE: Stepping stones are crucial for species’ long-distance dispersal and range expansion through habitat networks. Journal of Applied Ecology, 51(1), 171-182.\cr
#' -   Hanski, I. and Ovaskainen, O. 2000. The metapopulation capacity of a fragmented landscape. Nature 404: 755–758.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' data("habitat_nodes", package = "Makurhini")
#' nrow(habitat_nodes) # Number of patches
#' #Two distance threshold,
#' IIC <- MK_dPCIIC(nodes = habitat_nodes,
#'                 attribute = NULL,
#'                 area_unit = "m2",
#'                 distance = list(type = "centroid"),
#'                 LA = NULL,
#'                 overall = TRUE,
#'                 metric = "IIC",
#'                 distance_thresholds = c(10000, 20000)) #10,20 km
#' IIC$d20000$overall_d20000
#' plot(IIC$d20000$node_importances_d20000["dIIC"], breaks = "jenks")
#' plot(IIC$d20000$node_importances_d20000["dIICintra"], breaks = "jenks")
#' plot(IIC$d20000$node_importances_d20000["dIICflux"], breaks = "jenks")
#' plot(IIC$d20000$node_importances_d20000["dIICconnector"], breaks = "jenks")
#'
#' #Using raster and resistance
#' data("habitat_nodes_raster", package = "Makurhini")
#' data("resistance_matrix", package = "Makurhini")
#' PC <- MK_dPCIIC(nodes = habitat_nodes_raster,
#'                 attribute = NULL,
#'                 distance = list(type = "least-cost",
#'                                 resistance = resistance_matrix),
#'                 metric = "PC", probability = 0.5,
#'                 overall = TRUE,
#'                 distance_thresholds = 40000) # 40 km
#' PC$overall_d40000
#' PC$node_importances_d40000
#' plot(PC$node_importances_d40000)
#' }
#' @note Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @returns
#' -   If only \bold{one distance} was used in the parameter \code{distance_thresholds} then return an object of class \code{sf} with the node importance values (delta IIC or PC).\cr
#' -   If you add \bold{\code{overall = TRUE}}, then a list containing the \code{sf} class object with the importance values of the nodes and a \code{data.frame} with the overall connectivity values will be returned.\cr
#' -   If you add \bold{\code{overall = TRUE}} or \bold{\code{onlyoverall = TRUE}}, a \code{data.frame} with overall connectivity metrics will be returned, including \code{PC} or \code{IIC} (as selected), their corresponding \code{num} and \code{EC} values, and the fractions \code{PC/IIC intra} (intrapatch connectivity), \code{PC/IIC direct} (interpatch connectivity via direct links only), and \code{PC/IIC step} (interpatch connectivity attributable to stepping-stone use).\cr
#' -   If you use the \bold{\code{restoration}} parameter then an extra column will be returned to the \code{sf} object with the node importance values, unless you use the \code{onlyrestor} argument (i.e., equal to \code{TRUE}) only the restoration metric is estimated.\cr
#' -   If you use \bold{multiple distance thresholds} (e.g, \code{distance_thresholds = c(1000, 5000, 80000)}), the resulting data should be returned in the form of a \code{list}, wherein each \code{list} item contains the resulting objects for each distance threshold.
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom purrr map_dbl
#' @importFrom raster as.matrix extent raster stack extent<- writeRaster reclassify crs crs<- unique
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map
#' @importFrom igraph delete_vertices distances
#' @importFrom sf write_sf st_as_sf
#' @importFrom stats median
#' @importFrom terra unique rast classify writeRaster
MK_dPCIIC <- function(nodes,
                      attribute  = NULL,
                      weighted = FALSE,
                      LA = NULL,
                      area_unit = "m2",
                      restoration = NULL,
                      onlyrestor = FALSE,
                      distance = list(type= "centroid", resistance = NULL),
                      metric = c("IIC", "PC"),
                      probability = NULL,
                      distance_thresholds = NULL,
                      threshold = NULL,
                      overall = FALSE,
                      onlyoverall = FALSE,
                      parallel = NULL,
                      parallel_mode = 1,
                      write = NULL,
                      id_sel = NULL,
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

      if(is.null(probability)){
        probability = 0.5; message("Note: probability = NULL, then a probability of 0.5 will be used")
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
    parallel_mode <- 0
  } else if (isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
    if(any(parallel_mode == 0 | is.null(parallel_mode))){parallel_mode <- 1}
  } else if((!is.null(parallel))){
    if(!is.numeric(parallel)){
      stop("if you use parallel argument then you need a numeric value")
    }
    if(any(parallel_mode == 0 | is.null(parallel_mode))){parallel_mode <- 1}
  } else {
    parallel <- NULL
    parallel_mode <- 0
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
  } else if (any(class(nodes)[1] == "RasterLayer" | class(nodes)[1] == "SpatRaster")){
    rcl <- class(nodes)[1]
    if( rcl == "RasterLayer"){
      nodes <- terra::rast(nodes)
    }
    if(length(terra::unique(nodes)[[1]]) < 2){
      stop("error, you need more than 2 nodes")
    }
    idT <- NULL
    attribute_1 <- nodesfile(nodes, id = idT,
                             attribute = attribute,
                             area_unit = area_unit,
                             restoration = restoration,
                             weighted = weighted,
                             write = NULL)
    id_original <- attribute_1[,1]; attribute_1[,1] <- 1:nrow(attribute_1)
  } else if (class(nodes)[1] == "data.frame"){
    attribute_1 <- nodes; idT <- NULL; id_original <- attribute_1[,1]
    attribute_1[,1] <- 1:nrow(nodes); names(attribute_1)[1:2] <- c("IdTemp", "Area")

    if(dim(nodes)[2] < 2){
      stop("You need two or three columns, see the ?MK_dPCIIC")
    }

  } else {
    stop("error missing file of nodes")
  }

  if(!is.null(restoration)){
    if(isFALSE(ncol(attribute_1) == 3 & isTRUE(unique(unique(attribute_1[, 3]) <= 1)))){
      stop("error, review restoration argument, unique values must be 0 and 1")
    } else {
      names(attribute_1)[3] <- "res"
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
    rownames(dist) <- 1:nrow(dist); colnames(dist) <- 1:ncol(dist)
  } else if (class(distance)[1] == "data.frame"){
    if(dim(distance)[2] != 3){
      stop("If the distance is pairwise, it must have three columns: From, To, and Distance")
    }
  dist <- to_sym_matrix(distance,
                        from = names(distance)[1],
                        to = names(distance)[2],
                        value = names(distance)[3],
                        diag_value = 0, fill_na = 0)
  rownames(dist) <- 1:nrow(dist); colnames(dist) <- 1:ncol(dist)
  } else {
    if(isTRUE(intern)){
      if(!is.null(distance$resistance)){
        message("Estimating distances. This may take several minutes depending on the number of nodes and raster resolution")
      } else {
        if(nrow(attribute_1) > 1000){
          message("Estimating distances. This may take several minutes depending on the number of nodes")
        }
      }
    }
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

    if(is.null(idT)){
      rownames(dist) <- 1:nrow(dist); colnames(dist) <- 1:ncol(dist)
    }
  }

  if(is.null(distance_thresholds)){
    distance_thresholds <- median(dist)
  }

  if(isTRUE(intern)){
    message(paste0("Estimating ", metric, " index. This may take several minutes depending on the number of nodes"))
  }

  if(length(distance_thresholds) > 1 & isTRUE(intern)){
    pb <- txtProgressBar(0,length(distance_thresholds), style = 3)
  }
#x=1
  result <- lapply(1:length(distance_thresholds), function(x){
    x.1 <- distance_thresholds[x]

    if(isTRUE(threshold)){
      dist.i <- dist
      if(length(threshold) > 1){
        dist.i[which(dist.i <= threshold[x])] <- NA
      } else{
        dist.i[which(dist.i <= threshold)] <- NA
      }
      mat1 <- tryCatch(get_sdist(dist_nodes = dist.i,
                                 attr_nodes = attribute_1[,2],
                                 metric = metric,
                                 probability = probability,
                                 distance_threshold = x.1,
                                 igraph_Dijkstra = instr_dist,
                                 parallel = parallel,
                                 return_graph = TRUE,
                                 min_nodes = 0,
                                 loop = TRUE, G1 = 1000,
                                 pij_min = 0.01,
                                 return_pij = TRUE,
                                 intern = if(!is.null(parallel)){intern} else {FALSE}), error = function(err)err)
    } else {
      mat1 <- tryCatch(get_sdist(dist_nodes = dist,
                                 attr_nodes = attribute_1[,2],
                                 metric = metric,
                                 probability = probability,
                                 distance_threshold = x.1,
                                 igraph_Dijkstra = instr_dist,
                                 parallel = parallel,
                                 return_graph = TRUE,
                                 min_nodes = 0,
                                 loop = TRUE, G1 = 1000,
                                 pij_min = 0.01,
                                 return_pij = TRUE,
                                 intern = if(!is.null(parallel)){intern} else {FALSE}), error = function(err)err)
    }

    if(inherits(mat1, "error")){
      stop("Error in short distance estimation")
    } else {
      attribute_2 <- attribute_1[,2]
      num <- sum(mat1[[1]], na.rm = TRUE)
    }

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

      matr <- outer(attribute_2, attribute_2)
      intra <- sum(diag(matr)); diag(matr) <- 0

      if(metric == "PC"){
        direct <- sum(matr * mat1$pij, na.rm = TRUE)
        Step <- sum(matr * (mat1$`pij*`- mat1$pij), na.rm = TRUE)
      } else {
        matr2 <- mat1$`pij*`
        matr2[matr2 == 16443701] = 0
        direct <- sum(matr[which(matr2 == 1)]/2)
        Step <- sum(matr[which(matr2 >= 2)]/(1+ mat1$`pij*`[which(matr2 >= 2)]))
      }

      intra <- round((intra/num)*100, 5); direct <- round((direct/num)*100, 5)
      Step <- round((Step/num)*100, 5); ctrl <- sum(intra, direct, Step)

      if(ctrl < 100){
        ctrl <- (100 - ctrl)/3; intra <- intra + ctrl
        direct <- direct + ctrl; Step <- Step + ctrl
      }

      overall.2 <- rbind(overall.2,
                         data.frame(Index = paste0(metric, c("intra(%)", "direct(%)",
                                              "step(%)")),
                                    Value = c(intra, direct, Step)))
      if (!is.null(write)){
        write.csv(overall.2, file = paste0(write, "_Overall","_", "d", x.1,  ".csv"), row.names = FALSE)
      }
    }

    if(isFALSE(onlyoverall)){
      modo <- if(is.null(parallel_mode)){
        0
      } else {
        if(is.null(parallel)){0} else {parallel_mode}
      }

      if(length(attribute_2) > 500){
        if(isTRUE(intern)){
          message("Estimation of the index delta. Processing time depends on the number of patches.")
        }
        intern = TRUE
      }

      if(isFALSE(onlyrestor)){
        if(metric == "PC"){
          delta <- nodes_dPC(graph_n = mat1$graph,
                             attrib_n = attribute_2,
                             num = num,
                             G1 = 1000,
                             min_nodes = 0,
                             parallel = parallel,
                             parallel_mode = modo,
                             intern = intern)
        } else {
          delta <- nodes_dIIC(graph_n = mat1$graph,
                              attrib_n = attribute_2,
                              id_sel = id_sel,
                              num = num,
                              G1 = 1000,
                              min_nodes = 0,
                              parallel = parallel,
                              parallel_mode = modo,
                              intern = intern)
        }

        dintra <- round((((attribute_2^2)) / num) * 100, 7); dflux <- round(2*(rowSums(mat1[[1]], na.rm = TRUE) - attribute_2^2)/ num * 100, 7)
        dconnector <- round(map_dbl(delta - ((((attribute_2^2)) / num) * 100) - (2*(rowSums(mat1[[1]], na.rm = TRUE) - attribute_2^2)/ num * 100), function(y){if(y < 0){0} else {y}}), 20)

        metric_conn <- data.frame("IdTemp2" = if(is.null(idT)){id_original} else{attribute_1[,1]},
                                  "delta" = round(delta, 7),
                                  "dintra" = dintra, "dflux" = dflux, "dconnector" = dconnector,
                                  check.names = FALSE)
      } else {
        metric_conn <- data.frame("IdTemp2" = if(is.null(idT)){id_original} else{attribute_1[,1]},
                                  check.names = FALSE)
      }

      if(!is.null(restoration)){
        IdT0 <- attribute_1[which(attribute_1[,3] == 0), 1]
        if(metric == "PC"){
          graph_nodes.p1 <- mat1[[2]]
          graph_nodes.p1$data <- graph_nodes.p1$data[-which(graph_nodes.p1$data$from %in% (IdT0-1)|graph_nodes.p1$data$to %in% (IdT0-1)),]
          graph_nodes.p1$dict <- graph_nodes.p1$dict[-which(graph_nodes.p1$dict$ref %in% IdT0),]
          toV <- as.numeric(graph_nodes.p1$dict$ref)
          smat <- get_distance_matrix(Graph = graph_nodes.p1, from = toV, to = toV, algorithm = "phast")

          if(inherits(smat, "error")){
            stop("Error in short distance estimation")
          } else {
            smat <- exp(-smat); smat[is.infinite(smat)] <- 0
            smat <- outer(attribute_2[which(attribute_1[,3] == 1)], attribute_2[which(attribute_1[,3] == 1)]) * smat
            num_p1 <- sum(smat, na.rm = TRUE)
          }
          rm(graph_nodes.p1, smat, toV)

          dres <- rep(0, nrow(attribute_1))
          dres[IdT0] <- nodes_dPC_res(idres = IdT0,
                                      graph_n = mat1[[2]],
                                      attrib_n = attribute_2,
                                      num = num_p1,
                                      G1 = 1000,
                                      min_nodes = 0,
                                      parallel = parallel,
                                      parallel_mode = modo,
                                      intern = intern)
        } else {
          graph_nodes.p1 <- mat1[[2]]
          graph_nodes.p1 <- igraph::delete_vertices(graph_nodes.p1, IdT0)

          smat <- tryCatch(igraph::distances(graph_nodes.p1, weights = NULL),
                           error = function(err) err)

          if(inherits(smat, "error")){
            stop("Error in short distance estimation")
          } else {
            smat2 <- outer(attribute_2[which(attribute_1[,3] == 1)], attribute_2[which(attribute_1[,3] == 1)]) / (1 + smat); smat2[which(is.infinite(smat))] <- 0
            num_p1 <- sum(smat2, na.rm = TRUE)
          }
          rm(graph_nodes.p1, smat)

          dres <- rep(0, nrow(attribute_1))
          dres[IdT0] <- nodes_dIIC_res(idres = IdT0,
                                      graph_n = mat1[[2]],
                                      attrib_n = attribute_2,
                                      num = num_p1,
                                      G1 = 1000,
                                      min_nodes = 0,
                                      parallel = parallel,
                                      parallel_mode = modo,
                                      intern = intern)
          }

        metric_conn <- cbind(metric_conn, dres)
      } else {
        if(!is.null(parallel)){
          close_multiprocess(works)
        }
      }

      if(!is.null(idT)){
        nodes.2 <- cbind(nodes, metric_conn)
        if(isFALSE(onlyrestor)){
          names(nodes.2)[which(names(nodes.2) %in%
                                 c("delta", "dintra", "dflux", "dconnector"))] <- c(paste0("d", metric),
                                                                                    paste0("d", metric, "intra"),
                                                                                    paste0("d", metric, "flux"),
                                                                                    paste0("d", metric, "connector"))
        }

        if(!is.null(restoration)){
          names(nodes.2)[which(names(nodes.2) == "dres")] <- paste0("d", metric, "res")
        }

        nodes.2$IdTemp <- NULL; nodes.2$IdTemp2 <- NULL; nodes.2 <- nodes.2[,which(names(nodes.2) != "geometry")]

        if(!is.null(id_sel)){
          nodes.2 <- nodes.2[id_sel,]
        }

        if(!is.null(write)){
          write_sf(nodes.2, paste0(write, "_", "d", x.1,  ".shp"), delete_layer = TRUE)
        }

      } else {
        if(isFALSE(onlyrestor)){
          if(is.null(restoration)){
            names(metric_conn)[2:ncol(metric_conn)] <- c(paste0("d", metric),
                                                         paste0("d", metric, "intra"),
                                                         paste0("d", metric, "flux"),
                                                         paste0("d", metric, "connector"))
          } else {
            names(metric_conn)[2:ncol(metric_conn)] <- c(paste0("d", metric),
                                                         paste0("d", metric, "intra"),
                                                         paste0("d", metric, "flux"),
                                                         paste0("d", metric, "connector"),
                                                         paste0("d", metric, "res"))
          }
        } else {
          names(metric_conn)[2:ncol(metric_conn)] <- c(paste0("d", metric),
                                                       paste0("d", metric, "res"))
        }

        if(class(nodes)[1] == "data.frame"){
          nodes.2 <- cbind(nodes, metric_conn); nodes.2$IdTemp <- NULL; nodes.2$IdTemp2 <- NULL
          if(!is.null(id_sel)){
            nodes.2 <- nodes.2[id_sel,]
          }

          if(!is.null(write)){
            write.csv(nodes.2, paste0(write, "_", "d", x.1,  ".csv"), row.names = FALSE)
          }
        } else {
          rp <- unique(nodes)[[1]]; rp <- as.vector(rp); rp <- rp[which(!is.na(rp))]

          if(!is.null(parallel)){
            if(class(nodes)[1] == "SpatRaster"){
              nodes <- raster(nodes)
            }
            m <- matrix(nrow = nrow(attribute_1), ncol = 2); m[,1] <- id_original

            works <- as.numeric(availableCores())-1; works <-  if(parallel > works){works}else{parallel}
            if(.Platform$OS.type == "unix") {
              strat <- future::multicore
            } else {
              strat <- future::multisession
            }
            plan(strategy = strat, gc = TRUE, workers = works)

            r_metric <- tryCatch(future_map(2:ncol(metric_conn), function(c){
              x1 <- metric_conn[,c(1, c)]
              for(i in rp){
                m[which(m[,1] == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)}, .progress = TRUE), error = function(err) err)
            close_multiprocess(works)
            r_metric <- lapply(r_metric, rast); nodes <- rast(nodes)
          } else {
            m <- matrix(nrow = nrow(attribute_1), ncol = 2); m[,1] <- id_original
            r_metric <- lapply(2:ncol(metric_conn), function(c){
              x1 <- metric_conn[,c(1, c)]
              for(i in rp){
                m[which(m[,1] == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- classify(nodes, rcl = m)
              return(x1)})
          }

          nodes.2 <- list()
          for(i in 2:(length(r_metric)+1)){
            nodes.2[[i]] <- r_metric[[i-1]]
          }

          nodes.2[[1]] <- nodes; nodes.2 <- rast(nodes.2)
          names(nodes.2) <- names(metric_conn); names(nodes.2)[1] <- "Id"

          if (!is.null(write)){
            n <- names(metric_conn); n <- lapply(as.list(2:length(n)), function(w){
              x1 <- nodes.2[[w]]; crs(x1) <- crs(nodes.2)
              writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x.1, ".tif"),
                          overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
            })
          }

          if(rcl != "SpatRaster"){
            nodes.2 <- lapply(nodes.2, raster) |> stack(x = _)
          }
        }
      }

      if(isTRUE(overall)){
        result_metric <- list(nodes.2, overall.2)
        names(result_metric) <- c(paste0("node_importances_d",x.1), paste0("overall_d", x.1))
      } else {
        if(!is.null(id_sel)){
          metric_conn <- metric_conn[id_sel,]
        }
        result_metric <- nodes.2
      }

    } else {
      result_metric <- overall.2
    }

    if(length(distance_thresholds) > 1 & isTRUE(intern)){
      setTxtProgressBar(pb, x)
    }

    return(result_metric) })
  if(isTRUE(intern)){message("")}
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
  if(isTRUE(intern)){
    message("Done!")
  }
  return(result)
}

#' Estimate restoration using PC
#'
#' @param idres integer. Id patches = 0
#' @param graph_n Graph from cppRouting package
#' @param attrib_n numeric. Nodes attribute
#' @param num numeric. numPC ID patches = 1
#' @param G1 numeric
#' @param min_nodes numeric
#' @param parallel numeric
#' @param parallel_mode numeric
#' @param intern logical
#' @importFrom utils object.size
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map_dbl future_map_dfc
#' @importFrom purrr map_dbl map_dfc
#' @importFrom cppRouting get_distance_matrix
#' @keywords internal
nodes_dPC_res <- function(idres = NULL,
                          graph_n = NULL,
                          attrib_n = NULL,
                          num = NULL,
                          G1 = 1000,
                          min_nodes = 2000,
                          parallel = NULL,
                          parallel_mode = 0,
                          intern = FALSE){
  if(length(idres) < min_nodes){
    parallel_mode = 0
  }
  attribute_p1 <- attrib_n[-idres]; attribute_p0 <- attrib_n[idres]
  if(parallel_mode == 0){
    delta <- map_dbl(1:length(idres), function(i){
      graph_nodes.i <- graph_n; attribute.i <- c(attribute_p1, attribute_p0[i])
      graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from %in% (idres[-i]-1)|graph_nodes.i$data$to %in% (idres[-i]-1)),]
      graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref %in% idres[-i]),]

      toV <- as.numeric(graph_nodes.i$dict$ref)
      smat <- get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast")

      if(inherits(smat, "error")){
        stop("Error in short distance estimation")
      } else {
        smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
        num.i <- sum(smat, na.rm = TRUE); dC <- (num.i - num) / num.i * 100; return(dC)
      }
    }, .progress = intern)
    return(delta)
  } else {
    if(parallel_mode == 1){
      works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize = m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- future_map_dbl(1:length(idres), function(i){
        graph_nodes.i <- graph_n; attribute.i <- c(attribute_p1, attribute_p0[i])
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from %in% (idres[-i]-1)|graph_nodes.i$data$to %in% (idres[-i]-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref %in% idres[-i]),]

        toV <- as.numeric(graph_nodes.i$dict$ref)
        smat <- get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast")

        if(inherits(smat, "error")){
          stop("Error in short distance estimation")
        } else {
          smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
          num.i <- sum(smat, na.rm = TRUE); dC <- (num.i - num) / num.i * 100; return(dC)
        }
      }, .progress = intern); close_multiprocess(works)
      return(delta)
    } else {
      works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize= m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- map_dbl(1:length(idres), function(i){
        graph_nodes.i <- graph_n; attribute.i <- c(attribute_p1, attribute_p0[i])
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from %in% (idres[-i]-1)|graph_nodes.i$data$to %in% (idres[-i]-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref %in% idres[-i]),]

        toV <- as.numeric(graph_nodes.i$dict$ref)

        if(length(toV) >= min_nodes){
          seq_n <- seq(1,length(toV), G1)
          if(is.null(parallel)){
            smat <- tryCatch(map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- get_distance_matrix(Graph = graph_nodes.i,
                                                                    from = toV,
                                                                    to = toV.i[which(!is.na(toV.i))],
                                                                    allcores = FALSE,
                                                                    algorithm = "phast")
              return(y.1)
            }, .progress = intern)|> as.matrix(x = _), error = function(err)err)
          } else {
            smat <- tryCatch(furrr::future_map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- get_distance_matrix(Graph = graph_nodes.i,
                                                                    from = toV, to = toV.i[which(!is.na(toV.i))],
                                                                    allcores=FALSE, algorithm = "phast")
              invisible(gc())
              return(y.1)
            }, .progress = intern)|> as.matrix(x = _), error = function(err)err); close_multiprocess(works)

            if(inherits(smat, "error")){
              close_multiprocess(works)
              message("error probably due to the lack of memory in the parallel process. Makurhini will try using sequential process or you can change allocated memory (e.g., options(future.globals.maxSize= 2118123520)),
              before run again this function")
              smat <- tryCatch(map_dfc(seq_n, function(y){
                toV.i <- toV[y:(y + G1-1)];y.1 <- get_distance_matrix(Graph = graph_nodes.i,
                                                                      from = toV,
                                                                      to = toV.i[which(!is.na(toV.i))],
                                                                      allcores = FALSE,
                                                                      algorithm = "phast")
                return(y.1)
              }, .progress = intern)|> as.matrix(x = _), error = function(err)err)
            }
          }
        } else {
          smat <- tryCatch(get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast"), error = function(err)err)
        }

        if(inherits(smat, "error")){
          stop("Error in short distance estimation")
        } else {
          smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
          num.i <- sum(smat, na.rm = TRUE); dC <- (num.i - num) / num.i * 100; return(dC)
        }
      }, .progress = intern); close_multiprocess(works)
    }
  }
}

#' Estimate restoration using IIC
#'
#' @param idres integer. Id patches = 0
#' @param graph_n Graph from igraph package
#' @param attrib_n numeric. Nodes attribute
#' @param num numeric. numPC ID patches = 1
#' @param G1 numeric
#' @param min_nodes numeric
#' @param parallel numeric
#' @param parallel_mode numeric
#' @param intern logical
#' @importFrom utils object.size
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map_dbl future_map_dfc
#' @importFrom purrr map_dbl map_dfc
#' @importFrom igraph distances as_ids V delete_vertices
#' @keywords internal
nodes_dIIC_res <- function(idres = NULL,
                           graph_n = NULL,
                           attrib_n = NULL,
                           num = NULL,
                           G1 = 1000,
                           min_nodes = 2000,
                           parallel = NULL,
                           parallel_mode = 0,
                           intern = FALSE){
  if(length(idres) < min_nodes){
    parallel_mode = 0
  }
  attribute_p1 <- attrib_n[-idres]; attribute_p0 <- attrib_n[idres]
  if(parallel_mode == 0){
    delta <- map_dbl(1:length(idres), function(i) {
      graph_nodes.i <- graph_n; attribute.i <- c(attribute_p1, attribute_p0[i])

      delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, idres[-i]),
                          error = function(err) err)

      if(!inherits(delet.i, "error")){
        graph_nodes.i <- delet.i
      }

      smat.i <- tryCatch(igraph::distances(graph_nodes.i, weights = NULL),
                         error = function(err) err)

      if (inherits(smat.i, "error")) {
        stop("Error in short distance estimation")
      } else {
        smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
        num.i <- sum(smat2, na.rm = TRUE); dC <- (num.i - num)/num.i * 100
        return(dC)
      }
    }, .progress = intern)
    return(delta)
  } else {
    if(parallel_mode == 1){
      works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize = m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- future_map_dbl(1:length(idres), function(i) {
        graph_nodes.i <- graph_n; attribute.i <- c(attribute_p1, attribute_p0[i])

        delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, idres[-i]),
                            error = function(err) err)

        if(!inherits(delet.i, "error")){
          graph_nodes.i <- delet.i
        }

        smat.i <- tryCatch(igraph::distances(graph_nodes.i, weights = NULL),
                           error = function(err) err)

        if (inherits(smat.i, "error")) {
          stop("Error in short distance estimation")
        } else {
          smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
          num.i <- sum(smat2, na.rm = TRUE); dC <- (num.i - num)/num.i * 100
          return(dC)
        }
      }, .progress = intern); close_multiprocess(works)
      return(delta)
    } else {
      works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize= m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- map_dbl(1:length(idres), function(i) {
        graph_nodes.i <- graph_n; attribute.i <- c(attribute_p1, attribute_p0[i])

        delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, idres[-i]),
                            error = function(err) err)

        if(!inherits(delet.i, "error")){
          graph_nodes.i <- delet.i
        }
        toV <- as_ids(V(graph_nodes.i))

        if(length(toV) > min_nodes){
          seq_n <- seq(1,length(toV), G1)
          if(!is.null(parallel)){
            smat.i <- tryCatch(furrr::future_map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                  to = toV.i[which(!is.na(toV.i))],
                                                                  weights = NULL)
              return(y.1)}, .progress = intern)|> as.matrix(x = _),
              error = function(err)err)
            if(inherits(smat.i, "error")){
              close_multiprocess(works)
              message("error probably due to the lack of memory in the parallel process. Makurhini will try using sequential process or you can change allocated memory (e.g., options(future.globals.maxSize= 2118123520)), before run again this function")
              smat.i <- tryCatch(map_dfc(seq_n, function(y){
                toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                    to = toV.i[which(!is.na(toV.i))],
                                                                    weights = NULL)
                return(y.1)}, .progress = intern) |> as.matrix(x = _),
                error = function(err)err)
            }
          } else {
            smat.i <- tryCatch(map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                  to = toV.i[which(!is.na(toV.i))],
                                                                  weights = NULL)
              return(y.1)}, .progress = intern)|> as.matrix(x = _),
              error = function(err)err)
          }
        } else {
          smat.i <- tryCatch(igraph::distances(graph_nodes.i,
                                               weights = NULL),
                             error = function(err) err)
        }

        if (inherits(smat.i, "error")) {
          stop("Error in short distance estimation")
        } else {
          smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
          num.i <- sum(smat2, na.rm = TRUE); dC <- (num.i - num)/num.i * 100
          return(dC)
        }
      }, .progress = intern); close_multiprocess(works)
    }
  }
}

#' Estimate dPC
#'
#' @param graph_n Graph from cppRouting package
#' @param attrib_n numeric. Nodes attribute
#' @param num numeric. numPC ID patches = 1
#' @param G1 numeric
#' @param min_nodes numeric
#' @param parallel numeric
#' @param parallel_mode numeric
#' @param intern logical
#' @importFrom utils object.size
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map_dbl future_map_dfc
#' @importFrom purrr map_dbl map_dfc
#' @importFrom cppRouting get_distance_matrix
#' @keywords internal
nodes_dPC <- function(graph_n = NULL,
                      id_sel = NULL,
                      attrib_n = NULL,
                      num = NULL,
                      G1 = 1000,
                      min_nodes = 2000,
                      parallel = NULL,
                      parallel_mode = 0,
                      intern = FALSE){
  if(any(length(attrib_n) < min_nodes | !is.null(id_sel))){
    parallel_mode = 0
  }
  if(parallel_mode == 0){
    if(is.null(id_sel)){
      delta <- map_dbl(1:length(attrib_n), function(i){
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from == (i-1)|graph_nodes.i$data$to == (i-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref == i),]
        toV <- as.numeric(graph_nodes.i$dict$ref)
        smat <- get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast")

        if(inherits(smat, "error")){
          stop("Error in short distance estimation")
        } else {
          smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
          num.i <- sum(smat, na.rm = TRUE); dC <- (num - num.i) / num * 100; return(dC)
        }
      }, .progress = intern)
    } else {
      delta <- map_dbl(id_sel, function(i){
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from == (i-1)|graph_nodes.i$data$to == (i-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref == i),]
        toV <- as.numeric(graph_nodes.i$dict$ref)
        smat <- get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast")

        if(inherits(smat, "error")){
          stop("Error in short distance estimation")
        } else {
          smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
          num.i <- sum(smat, na.rm = TRUE); dC <- (num - num.i) / num * 100; return(dC)
        }
      }, .progress = intern)
    }

    return(delta)
  } else {
    if(parallel_mode == 1){
      works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize = m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- future_map_dbl(1:length(attrib_n), function(i){
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from == (i-1)|graph_nodes.i$data$to == (i-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref == i),]
        toV <- as.numeric(graph_nodes.i$dict$ref)
        smat <- get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast")

        if(inherits(smat, "error")){
          stop("Error in short distance estimation")
        } else {
          smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
          num.i <- sum(smat, na.rm = TRUE); dC <- (num - num.i) / num * 100; return(dC)
        }
       }, .progress = intern); close_multiprocess(works)
      return(delta)
    } else {
      works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize= m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- map_dbl(1:length(attrib_n), function(i){
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        graph_nodes.i$data <- graph_nodes.i$data[-which(graph_nodes.i$data$from == (i-1)|graph_nodes.i$data$to == (i-1)),]
        graph_nodes.i$dict <- graph_nodes.i$dict[-which(graph_nodes.i$dict$ref == i),]
        toV <- as.numeric(graph_nodes.i$dict$ref)

        if(length(toV) >= min_nodes){
          seq_n <- seq(1,length(toV), G1)
          if(is.null(parallel)){
            smat <- tryCatch(map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- get_distance_matrix(Graph = graph_nodes.i,
                                                                    from = toV,
                                                                    to = toV.i[which(!is.na(toV.i))],
                                                                    allcores = FALSE,
                                                                    algorithm = "phast")
              return(y.1)
            }, .progress = intern)|> as.matrix(x = _), error = function(err)err)
          } else {
            smat <- tryCatch(furrr::future_map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- get_distance_matrix(Graph = graph_nodes.i,
                                                                    from = toV, to = toV.i[which(!is.na(toV.i))],
                                                                    allcores=FALSE, algorithm = "phast")
              invisible(gc())
              return(y.1)
            }, .progress = intern)|> as.matrix(x = _), error = function(err)err); close_multiprocess(works)

            if(inherits(smat, "error")){
              close_multiprocess(works)
              message("error probably due to the lack of memory in the parallel process. Makurhini will try using sequential process or you can change allocated memory (e.g., options(future.globals.maxSize= 2118123520)),
              before run again this function")
              smat <- tryCatch(map_dfc(seq_n, function(y){
                toV.i <- toV[y:(y + G1-1)];y.1 <- get_distance_matrix(Graph = graph_nodes.i,
                                                                      from = toV, to = toV.i[which(!is.na(toV.i))],
                                                                      allcores=FALSE, algorithm = "phast")
                return(y.1)
              }, .progress = intern)|> as.matrix(x = _), error = function(err)err)
            }
          }
        } else {
          smat <- tryCatch(get_distance_matrix(Graph = graph_nodes.i, from = toV, to = toV, algorithm = "phast"), error = function(err)err)
        }

        if(inherits(smat, "error")){
          stop("Error in short distance estimation")
        } else {
          smat <- exp(-smat); smat[is.infinite(smat)] <- 0; smat <- outer(attribute.i, attribute.i) * smat
          num.i <- sum(smat, na.rm = TRUE); dC <- (num - num.i) / num * 100; return(dC)
        }
      }, .progress = intern); close_multiprocess(works)
    }
  }
}

#' Estimate dIIC
#'
#' @param graph_n Graph from igraph package
#' @param attrib_n numeric. Nodes attribute
#' @param num numeric. numPC ID patches = 1
#' @param G1 numeric
#' @param min_nodes numeric
#' @param parallel numeric
#' @param parallel_mode numeric
#' @param intern logical
#' @importFrom utils object.size
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map_dbl future_map_dfc
#' @importFrom purrr map_dbl map_dfc
#' @importFrom igraph distances as_ids V delete_vertices
#' @keywords internal
nodes_dIIC <- function(graph_n = NULL,
                       id_sel = NULL,
                       attrib_n = NULL,
                       num = NULL,
                       G1 = 1000,
                       min_nodes = 2000,
                       parallel = NULL,
                       parallel_mode = 0,
                       intern = FALSE){
  if(any(length(attrib_n) < min_nodes | !is.null(id_sel))){
    parallel_mode = 0
  }

  if(parallel_mode == 0){
    if(is.null(id_sel)){
      delta <- map_dbl(1:length(attrib_n), function(i) {
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, i),
                            error = function(err) err)

        if(!inherits(delet.i, "error")){
          graph_nodes.i <- delet.i
        }

        smat.i <- tryCatch(igraph::distances(graph_nodes.i,
                                             weights = NULL),
                           error = function(err) err)

        if (inherits(smat.i, "error")) {
          stop("Error in short distance estimation")
        } else {
          smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
          num.i <- sum(smat2, na.rm = TRUE); dC <- (num - num.i)/num * 100
          return(dC)
        }
      }, .progress = intern)
    } else {
     delta <-   map_dbl(id_sel, function(i) {
       graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
       delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, i),
                           error = function(err) err)

       if(!inherits(delet.i, "error")){
         graph_nodes.i <- delet.i
       }

       smat.i <- tryCatch(igraph::distances(graph_nodes.i,
                                            weights = NULL),
                          error = function(err) err)

       if (inherits(smat.i, "error")) {
         stop("Error in short distance estimation")
       } else {
         smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
         num.i <- sum(smat2, na.rm = TRUE); dC <- (num - num.i)/num * 100
         return(dC)
       }
     }, .progress = intern)
    }
    return(delta)
  } else {
    if(parallel_mode == 1){
      works <- as.numeric(availableCores())-1
      works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize = m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- future_map_dbl(1:length(attrib_n), function(i) {
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, i),
                            error = function(err) err)

        if(!inherits(delet.i, "error")){
          graph_nodes.i <- delet.i
        }
        toV <- as_ids(V(graph_nodes.i))

        if(length(toV) > min_nodes){
          seq_n <- seq(1,length(toV), G1)
          smat.i <- tryCatch(map_dfc(seq_n, function(y){
            toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                to = toV.i[which(!is.na(toV.i))],
                                                                weights = NULL)
            return(y.1)}, .progress = intern)|> as.matrix(x = _), error = function(err)err)

        } else {
          smat.i <- tryCatch(igraph::distances(graph_nodes.i,
                                               weights = NULL),
                             error = function(err) err)
        }

        if (inherits(smat.i, "error")) {
          stop("Error in short distance estimation")
        } else {
          smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
          num.i <- sum(smat2, na.rm = TRUE); dC <- (num - num.i)/num * 100
          return(dC)
        }
      }, .progress = intern); close_multiprocess(works)

      return(delta)
    } else {
      works <- as.numeric(availableCores())-1;works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      m <- as.numeric(object.size(graph_n))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize= m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
      delta <- map_dbl(1:length(attrib_n), function(i) {
        graph_nodes.i <- graph_n; attribute.i <- attrib_n[-i]
        delet.i <- tryCatch(igraph::delete_vertices(graph_nodes.i, i),
                            error = function(err) err)

        if(!inherits(delet.i, "error")){
          graph_nodes.i <- delet.i
        }
        toV <- as_ids(V(graph_nodes.i))

        if(length(toV) > min_nodes){
          seq_n <- seq(1,length(toV), G1)
          if(!is.null(parallel)){
            smat.i <- tryCatch(furrr::future_map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                  to = toV.i[which(!is.na(toV.i))],
                                                                  weights = NULL)
              return(y.1)}, .progress = intern)|> as.matrix(x = _),
              error = function(err)err)
            if(inherits(smat.i, "error")){
              close_multiprocess(works)
              message("error probably due to the lack of memory in the parallel process. Makurhini will try using sequential process or you can change allocated memory (e.g., options(future.globals.maxSize= 2118123520)), before run again this function")
              smat.i <- tryCatch(map_dfc(seq_n, function(y){
                toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                    to = toV.i[which(!is.na(toV.i))],
                                                                    weights = NULL)
                return(y.1)}, .progress = intern) |> as.matrix(x = _), error = function(err)err)
            }
          } else {
            smat.i <- tryCatch(map_dfc(seq_n, function(y){
              toV.i <- toV[y:(y + G1-1)];y.1 <- igraph::distances(graph_nodes.i,
                                                                  to = toV.i[which(!is.na(toV.i))],
                                                                  weights = NULL)
              return(y.1)}, .progress = intern)|> as.matrix(x = _), error = function(err)err)
          }
        } else {
          smat.i <- tryCatch(igraph::distances(graph_nodes.i,
                                               weights = NULL),
                             error = function(err) err)
        }

        if (inherits(smat.i, "error")) {
          stop("Error in short distance estimation")
        } else {
          smat2 <- outer(attribute.i, attribute.i) / (1 + smat.i); smat2[which(is.infinite(smat.i))] <- 0
          num.i <- sum(smat2, na.rm = TRUE); dC <- (num - num.i)/num * 100
          return(dC)
        }
      }, .progress = intern); close_multiprocess(works)
      return(delta)
    }
  }
}
