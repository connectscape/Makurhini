#' @title Estimate the probability of connectivity (PC) and prioritize nodes using Spatial Absorbing Markov Chains
#'
#' @description
#' Computes overall landscape connectivity and the importance (contribution) of
#' each node (habitat patch) to that connectivity using the Probability of
#' Connectivity (PC) index for a specified dispersal distance. Movement and
#' settlement are modeled under a Spatial Absorbing Markov Chain (SAMC)
#' framework.
#'
#' @details
#' The SAMC framework treats raster cells as states with three roles: departure, transient
#' (movement) and absorbing (settlement at patch boundaries). Transition
#' probabilities are derived from a resistance surface (and a chosen neighborhood
#' structure), while absorption rates reflect species’ dispersal and patch
#' properties. Multi-step transitions yield probabilities of reaching (and
#' settling in) patches across the landscape. These probabilities are then used
#' to compute the \strong{overall connectivity} (PC) and the \strong{node
#' importance} (each patch’s contribution to PC, e.g., via ΔPC). The SAMC
#' computations are implemented via the \pkg{samc} R package (Marx et al., 2020).
#' @param nodes Object containing nodes (e.g., habitat patches or fragments) information. It can be of the following classes:\cr
#' -   \code{Data.frame} with at least two columns: the first for node IDs and the second for attributes. If the `restoration` argument is used, the data frame must include a third column for restoration values.\cr
#' -   Spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.\cr
#' -   Raster (class \code{RasterLayer, SpatRaster}). It must be in a projected coordinate system. The values must be integers representing the ID of each habitat patch or node, with non-habitat areas represented by NA values (see \link[raster]{clump} or \link[terra]{patches}).
#' @param keep numeric. Use this option to simplify the nodes geometry and reduce the number of vertices. The value can range from 0 to 1 and is the proportion of points to retain (default equal to NULL). The higher the value, the higher the speed but less precision.
#' @param attribute \code{Character} or \code{vector}. If \code{NULL} (only applicable when \code{nodes} is of spatial data of vector or raster type) the area of the nodes will be used as the node attribute. The unit of area can be selected using the \code{area_unit} parameter. To use an alternative attribute, consider the class type of the object in the \code{nodes} parameter: \cr
#' -   If \code{nodes} is a spatial vector or data.frame, specify \bold{the name of the column} containing the attribute for the nodes. \cr
#' -   If \code{nodes} is a raster layer then it must be a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes. If the parameter \bold{weighted} is \code{TRUE} then the numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#' @param weighted \code{logical}. If the \code{nodes} are of raster type, you can weight the estimated area of each node by the attribute. When using this parameter the \code{attribute} parameter, which must be a vector of length equal to the number of nodes, usually has values between 0 and 1.
#' @param LA \code{numeric}. (\emph{optional, default = } \code{NULL}). The maximum landscape attribute, which is the attribute value that would correspond to a hypothetical habitat patch covering all the landscape with the best possible habitat, in which PC would be equal to 1. For example, if nodes attribute corresponds to the node area, then LA equals total landscape area. If nodes attribute correspond to a quality-weighted area and the quality factor ranges from 0 to 100, LA will be equal to 100 multiplied by total landscape area. The value of LA does not affect at all the importance of the nodes and is only used to calculate the overall landscape connectivity. If no LA value is entered (default) and  \code{overall = TRUE} or \code{onlyoverall = TRUE}, the function will only calculate the numerator of the global connectivity indices and the equivalent connected ECA or EC index.
#' @param area_unit \code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.
#' @param restoration \code{Character or vector}, (\emph{optional}). This parameter specifies the binary restoration value indicating whether each node is existing or hypothetically added for restoration: 1 for existing nodes in the landscape and 0 for nodes to be added to the initial landscape. If \code{NULL} (default), all nodes are considered existing (as if all restoration values were 1). Otherwise: \cr
#' -   If \code{nodes} is a shapefile (spatial vector) or a data.frame, specify the name of the column containing the restoration values. If `nodes` is a data.frame, it must have three columns: the first for node IDs, the second for attributes, and the third for restoration values.\cr
#' -   If \code{nodes} is a raster layer, provide a numeric vector with restoration values for each node in the raster.
#' @param onlyrestor \code{logical}. If \code{TRUE}, then only restoration metric will be calculated.
#' @param disp_kernel
#' A \code{list} with the parameters used to estimate settlement probability.
#' Expected components:
#' \describe{
#'   \item{\code{resistance}}{Resistance raster (e.g., \code{terra::SpatRaster})
#'   where higher values indicate greater movement cost. The resistance raster must be a closed matrix with no isolated pixels.}
#'
#'   \item{\code{CostFun}}{Function to compute the cost of moving between cells
#'   from local resistance. Default: \code{function(x) 1/mean(x)}.}
#'
#'   \item{\code{absorption}}{raster with absorption rate for absorbing states. By default,
#'   it is estimated from the species’ mean dispersal distance and the cell size
#'   of the resistance raster (Fletcher et al., 2023). If raster, then it must be a closed matrix with no isolated pixels.}
#'
#'   \item{\code{dir}}{Neighborhood (directions) for the probability graph.
#'   Options: \code{4}, \code{8} (Moore, default), or \code{16}.}
#'
#'   \item{\code{SettlFunc}}{Function that aggregates (combines) settlement
#'   probabilities among polygon vertices (edges) for each patch. Default:
#'   \code{mean}.}
#'
#'   \item{\code{pvertices}}{Proportion of boundary vertices per patch to include
#'   in the analysis (range \code{(0, 1]}). Default: \code{0.1}.}
#'
#'    \item{\code{nvertices}}{Number of boundary vertices per patch to include
#'   in the analysis. Default: \code{0.1}.}
#'
#'    \item{\code{systematic}}{If TRUE, vertices are selected systematically, with each selected vertex separated by approximately the same number of intermediate vertices. When FALSE, the function selects vertices at random. Default: \code{TRUE}.}
#' }
#' @param probability A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{dispersal_dist}. For example, if the \code{dispersal_dist} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{dispersal_dist} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the \code{"PC"} metric. If \code{probability = NULL}, then a probability of 0.5 will be used.
#' @param dispersal_dist A \code{numeric} indicating the dispersal distance (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range.
#' @param overall \code{logical}. If \code{TRUE}, then the overall metrics will be added to the result which is transformed into a list. Default equal to FALSE
#' @param onlyoverall \code{logical}. If \code{TRUE}, then only overall metrics will be calculated.
#' @param parallel  (\emph{optional, default =} \code{NULL}).
#' A \code{numeric} specifying the number of cores to parallelize the index estimation of the PC index and its deltas.Particularly useful when you have more than 1000 nodes. By default the analyses are not parallelized.
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
#' -   Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.\cr
#' -   Saura, S., Bodin, Ö., & Fortin, M.-J. (2014). EDITOR’S CHOICE: Stepping stones are crucial for species’ long-distance dispersal and range expansion through habitat networks. Journal of Applied Ecology, 51(1), 171-182.\cr
#' -   Hanski, I. and Ovaskainen, O. 2000. The metapopulation capacity of a fragmented landscape. Nature 404: 755–758.\cr
#' -   Marx, A. J., Wang, C., Sefair, J. A., Acevedo, M. A., & Fletcher Jr., R. J. (2020). samc: An R package for connectivity modeling with spatial absorbing Markov chains. Ecography, 43(4), 518-527.\cr
#' -   Fletcher, R. J., Iezzi, M. E., Guralnick, R., Marx, A. J., Ryan, S. J., & Valle, D. (2023). A framework for linking dispersal biology to connectivity across landscapes. Landscape Ecology, 38(10), 2487-2500. https://doi.org/10.1007/s10980-023-01741-8
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(sf)
#' library(raster)
#' data("habitat_nodes", package = "Makurhini")
#' nrow(habitat_nodes) # Number of patches
#' data("resistance_matrix", package = "Makurhini")
#'
#'
#' #Example
#' hn <- habitat_nodes[350:404,]
#' hn$Id <- 1:nrow(hn)
#' plot(hn)
#' mask_hn <- st_as_sfc(st_bbox(hn)) |> st_as_sf() |> st_buffer(5000)
#' #Aggregate for this example
#' resist <- aggregate(resistance_matrix, 4)
#' resist <- crop(resistance_matrix, mask_hn)
#' plot(resist)
#'
#' #Two distance threshold,
#' PCsamc <- MK_dPC_SAMC(nodes = hn,
#'                       keep = 0.5,
#'                       attribute = NULL,
#'                       area_unit = "m2",
#'                       disp_kernel = list(resistance = resist,
#'                                          CostFun = function(x) 1/mean(x),
#'                                          absorption = NULL, dir = 4,
#'                                          SettlFunc = sum, nvertices = 10),
#'                       LA = NULL,
#'                       overall = TRUE,
#'                       dispersal_dist = 10000) #10 km
#' PCsamc$overall_d10000
#' plot(PCsamc$node_importances_d10000["dPC"], breaks = "jenks")
#' plot(PCsamc$node_importances_d10000["dPCintra"], breaks = "jenks")
#' plot(PCsamc$node_importances_d10000["dPCflux"], breaks = "jenks")
#' plot(PCsamc$node_importances_d10000["dPCconnector"], breaks = "jenks")
#'
#' }
#' @note Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @returns
#' -   If you add \bold{\code{overall = TRUE}}, then a list containing the \code{sf} class object with the importance values of the nodes and a \code{data.frame} with the overall connectivity values (num, EC, PC) will be returned.\cr
#' -   If you use the \bold{\code{restoration}} parameter then an extra column will be returned to the \code{sf} object with the node importance values, unless you use the \code{onlyrestor} argument (i.e., equal to \code{TRUE}) only the restoration metric is estimated.\cr
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom purrr map_dbl
#' @importFrom raster as.matrix extent raster stack extent<- writeRaster reclassify crs crs<- unique
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map
#' @importFrom igraph delete_vertices distances
#' @importFrom sf write_sf st_as_sf
#' @importFrom stats median
#' @importFrom terra unique rast classify writeRaster res
#' @importFrom samc samc
#' @importFrom rmapshaper ms_simplify

MK_dPC_SAMC <- function(nodes,
                        keep = NULL,
                        attribute  = NULL,
                        weighted = FALSE,
                        LA = NULL,
                        area_unit = "m2",
                        restoration = NULL,
                        onlyrestor = FALSE,
                        disp_kernel = list(resistance = NULL, CostFun = NULL,
                                           absorption = NULL, dir = 8,
                                           SettlFunc = sum,
                                           pvertices = NULL,
                                           nvertices = 10,
                                           systematic = TRUE),
                        probability = NULL,
                        dispersal_dist = NULL,
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

  instr_dist <- FALSE

  if(is.null(dispersal_dist)){
    stop("dispersal_dist missing")
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

  #samc
  if(is.null(disp_kernel)){
    stop("Missing disp_kernel")
  } else {
    if(is.null(disp_kernel$resistance)){
      stop("Missing resistance")
    } else {
      if(class(disp_kernel$resistance)[1] == "RasterLayer"){
        disp_kernel$resistance <- rast(disp_kernel$resistance)
      }
    }

    if(is.null(disp_kernel$CostFun)){
      disp_kernel$CostFun <- function(x) 1/mean(x)
    }

    if(is.null(disp_kernel$absorption)){
      habitat_nodes$valTemp <- 0
      hab_rast <- terra::rasterize(x = vect(habitat_nodes), y = disp_kernel$resistance, field = "valTemp")
      hab_rast[is.na(hab_rast)] <- 1
      dist_pixel <- dispersal_dist/res(disp_kernel$resistance)[1]
      rs <- exp(-1.247 - 1.788 * log(dist_pixel))
      disp_kernel$absorption <- hab_rast * rs
    }

    if(is.null(disp_kernel$dir)){
      disp_kernel$dir <- 8
    }

    if(is.null(disp_kernel$SettlFunc)){
      disp_kernel$SettlFunc <- mean
    }

    if(is.null(disp_kernel$vertices)){
      disp_kernel$vertices <- 0.1
    }
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
      stop("You need two or three columns, see the ?MK_dPC_SAMC")
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

  #Dispersal kernel
  samc_stack <- c(disp_kernel$resistance, disp_kernel$absorption)
  samc_stack[any(is.na(samc_stack))] <- NA

  message("Computing settlement probabilities. Duration varies with raster resolution and node count.")
  #samc class
  samc_obj <- samc(data = samc_stack[[1]], absorption = samc_stack[[2]],
                   model = list(fun = disp_kernel$CostFun,
                                dir = disp_kernel$dir,
                                sym = TRUE))

  pij <- samc_settlement(hn = nodes, samc_class = samc_obj,
                         pvertices = disp_kernel$pvertices,
                         nvertices = disp_kernel$nvertices,
                         SettlFunc = disp_kernel$SettlFunc,
                         sampletype = disp_kernel$systematic,
                         simpl = keep,
                         intern = intern)

  if(!is.null(idT)){
    rownames(pij) <- nodes$IdTemp; colnames(pij) <- nodes$IdTemp
  }

  ###
  if(isTRUE(intern)){
    message("Estimating PC index. This may take several minutes depending on the number of nodes")
  }

  ###
  mat1 <- tryCatch(get_sdist(dist_nodes = NULL,
                             pij_mat = pij,
                             attr_nodes = attribute_1[,2],
                             metric = "PC",
                             probability = probability,
                             distance_threshold = dispersal_dist,
                             igraph_Dijkstra = FALSE,
                             parallel = parallel,
                             return_graph = TRUE,
                             min_nodes = 0,
                             loop = TRUE, G1 = 1000,
                             pij_min = 0.001,
                             return_pij = TRUE,
                             intern = if(!is.null(parallel)){intern} else {FALSE}), error = function(err)err)

  if(inherits(mat1, "error")){
    stop("Error in short distance estimation")
  } else {
    attribute_2 <- attribute_1[,2]
    num <- sum(mat1[[1]], na.rm = TRUE)
  }

  if(isTRUE(overall) | isTRUE(onlyoverall)){
    if(!is.null(LA)){
      if(LA > sum(attribute_2)){
        overall.2 <- data.frame(Index = c("PCnum", "EC(PC)", "PC"),
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
      overall.2 <- data.frame(Index = c("PCnum", "EC(PC)"),
                              Value = c(num, sqrt(num)))
    }

    if (!is.null(write)){
      write.csv(overall.2, file = paste0(write, "_Overall","_", "d", dispersal_dist,  ".csv"), row.names = FALSE)
    }
  }
  ##
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
      delta <- nodes_dPC(graph_n = mat1$graph,
                         attrib_n = attribute_2,
                         num = num,
                         G1 = 1000,
                         min_nodes = 0,
                         parallel = parallel,
                         parallel_mode = modo,
                         intern = intern)
      dintra <- round((((attribute_2^2)) / num) * 100, 7); dflux <- round(2*(rowSums(mat1[[1]], na.rm = TRUE) - attribute_2^2)/ num * 100, 7)
      dconnector <- round(map_dbl(delta - ((((attribute_2^2)) / num) * 100) - (2*(rowSums(mat1[[1]], na.rm = TRUE) - attribute_2^2)/ num * 100), function(y){if(y < 0){0} else {y}}), 20)

      metric_conn <- data.frame("IdTemp2" = if(is.null(idT)){id_original} else{attribute_1[,1]},
                                "dPC" = round(delta, 7),
                                "dPCintra" = dintra, "dPCflux" = dflux, "dPCconnector" = dconnector,
                                check.names = FALSE)
    } else {
      metric_conn <- data.frame("IdTemp2" = if(is.null(idT)){id_original} else{attribute_1[,1]},
                                check.names = FALSE)
    }

    if(!is.null(restoration)){
      IdT0 <- attribute_1[which(attribute_1[,3] == 0), 1]

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

      names(dres) <- "dPCres"
      metric_conn <- cbind(metric_conn, dres)
    } else {
      if(!is.null(parallel)){
        close_multiprocess(works)
      }
    }

    if(!is.null(idT)){
      nodes.2 <- cbind(nodes, metric_conn)
      nodes.2$IdTemp <- NULL; nodes.2$IdTemp2 <- NULL; nodes.2 <- nodes.2[,which(names(nodes.2) != "geometry")]

      if(!is.null(id_sel)){
        nodes.2 <- nodes.2[id_sel,]
      }

      if(!is.null(write)){
        write_sf(nodes.2, paste0(write, "_", "d", dispersal_dist,  ".shp"), delete_layer = TRUE)
      }

    } else {
      if(class(nodes)[1] == "data.frame"){
        nodes.2 <- cbind(nodes, metric_conn); nodes.2$IdTemp <- NULL; nodes.2$IdTemp2 <- NULL
        if(!is.null(id_sel)){
          nodes.2 <- nodes.2[id_sel,]
        }

        if(!is.null(write)){
          write.csv(nodes.2, paste0(write, "_", "d", dispersal_dist,  ".csv"), row.names = FALSE)
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
            writeRaster(x1, filename = paste0(write, "_", n[w], "_",  dispersal_dist, ".tif"),
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
      names(result_metric) <- c(paste0("node_importances_d", dispersal_dist),
                                paste0("overall_d", dispersal_dist))
    } else {
      if(!is.null(id_sel)){
        metric_conn <- metric_conn[id_sel,]
      }
      result_metric <- nodes.2
    }

  } else {
    result_metric <- overall.2
  }

  if(isTRUE(intern)){message("")}
    if(isTRUE(intern)){
    message("Done!")
  }
  return(result_metric)
}


#' Estimate probability of settlement
#'
#' @param hn nodes
#' @param samc_class samc object
#' @param vertices numeric
#' @param SettlFunc function
#' @param sampletype logical
#' @param simpl numeric
#' @param intern logical
#' @importFrom samc locate mortality pairwise
#' @importFrom purrr compact
#' @importFrom rmapshaper ms_simplify
#' @keywords internal
samc_settlement <- function(hn = NULL,
                            samc_class = NULL,
                            pvertices  = NULL,
                            nvertices  = NULL,
                            SettlFunc = NULL,
                            sampletype = TRUE,
                            simpl = NULL,
                            intern = TRUE
){
  if(!is.null(simpl)){
    hn <- ms_simplify(input = hn, keep = simpl, keep_shapes = TRUE, explode = FALSE)
  }
  origins_dest <- lapply(1:nrow(hn), function(x){
    x.1 <- SR_vertices(poly = hn[x,], k = nvertices, kp = pvertices, systematic = sampletype)
    x.2 <- tryCatch(lapply(1:nrow(x.1), function(i){
      i.1 <- tryCatch(locate(samc_class, x.1[i,1:2]), error = function(err)err)
      if(inherits(i.1, "error")){return(NULL)}
      return(i.1)
    }), error = function(err)err)

    if(inherits(x.2, "error")){
      x.2 <- tryCatch(locate(samc_class, x.1[1, 1:2]), error = function(err)err)
    }
    x.2 <- compact(x.2); x.2 <- unlist(x.2, use.names = FALSE)
    return(x.2)
  })
  n <- length(origins_dest)

  if(isTRUE(intern)){
    pb <- txtProgressBar(0,n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  ###Settlement prob
  disp_matrix <- matrix(0, nrow = nrow(hn), ncol = nrow(hn))

  for(x in 1:n){
    x.Orig <- origins_dest[[x]]
    x.Dest <- setdiff(seq_len(n), x)
    x.Dest2 <- unlist(origins_dest[x.Dest], use.names = FALSE)

    # skip if no valid origin or dest cells
    if (length(x.Orig) == 0L || length(x.Dest2) == 0L) next

    x.2 <- pairwise(mortality, samc_class, origin = x.Orig, dest = x.Dest2)
    x.1 <- vapply(seq_along(x.Dest), function(k) {
      idx_k <- origins_dest[[x.Dest[k]]]
      # rows of pw that correspond to this patch's cells
      rows_k <- x.2[, 2] %in% idx_k
      SettlFunc(x.2[rows_k, 3])
    }, numeric(1))

    disp_matrix[x, x.Dest] <- x.1
    disp_matrix[x.Dest, x] <- x.1 # symmetric
    if(isTRUE(intern)){
      setTxtProgressBar(pb, x)
    }
  }
  diag(disp_matrix) <- 1
  rownames(disp_matrix) <- 1:nrow(hn)
  colnames(disp_matrix) <- 1:nrow(hn)
  return(disp_matrix)
}


#' Get systematic or random vertices
#'
#' @param poly nodes
#' @param k numeric Number of vertices
#' @param kp numeric Proportion of vertices
#' @param systematic logical
#' @param return_sf logical
#' @importFrom sf st_coordinates st_geometry st_as_sf st_crs
#' @keywords internal
SR_vertices <- function(poly = NULL, k = NULL, kp = NULL, systematic = TRUE, return_sf = TRUE) {
  # get geometry (first feature for simplicity)
  g <- st_geometry(poly)[[1]]

  # get coordinates (we assume one exterior ring)
  coords <- st_coordinates(g)[, c("X", "Y")]

  # polygons repeat the first point at the end -> drop last
  if (nrow(coords) > 1 && all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- coords[-nrow(coords), , drop = FALSE]
    coords <- coords[!duplicated(coords),]
  }

  n <- nrow(coords)

  if(!is.null(kp)){
    k <- round(nrow(coords) * kp)
  }

  if (k >= n) {
    # if asking for more/equal than existing vertices, just return all
    idx <- seq_len(n)
  } else {
    if(isTRUE(systematic)){
      # systematic selection of indices along 1..n
      idx <- unique(round(seq(1, n, length.out = k)))
    } else {
      set.seed(234); idx <- sample(1:n, k)
    }
  }

  out <- data.frame(
    vertex_id = idx,
    x = coords[idx, 1],
    y = coords[idx, 2]
  )

  if (!return_sf){
    return(out)
  } else {
    return(st_as_sf(out,
                    coords = c("x", "y"),
                    crs = st_crs(poly)))
    }
}
