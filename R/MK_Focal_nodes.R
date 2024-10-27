#' Estimate the focal integral index of connectivity or the focal probability of connectivity
#'
#' This function enables the calculation of the focal Integral Index of Connectivity (\eqn{\mathbf{IIC}_f}) or the focal Probability of Connectivity (\eqn{\mathbf{PC}_f}) under one or more distance thresholds. Furthermore, this function estimates the composite connectivity index (\eqn{\mathbf{CCI}_f}; for further details, please see Latorre-Cárdenas et al., 2023).
#' @param nodes \code{sf, SpatVector, SpatialPolygonsDataFrame}. Object containing nodes (e.g., habitat nodes or fragments) to analyze information. Nodes are spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.
#' @param id  \code{character}. Column name with the node ID.
#' @param attribute \code{character} or \code{vector}. If \code{NULL} the area of the nodes will be used as the node attribute. The unit of area can be selected using the \code{area_unit} parameter. To use an alternative attribute, consider the class type of the object in the \code{nodes} parameter: \cr
#' -   If \code{nodes} is a spatial vector or data.frame, specify \bold{the name of the column} containing the attribute for the nodes. \cr
#' @param raster_attribute \code{raster, rast}. Raster object used to assign the attribute values of each node using the function specified in the \bold{fun_attribute} parameter.
#' @param fun_attribute \code{function}. Specifies the function to estimate the node attribute when the \bold{raster_attribute} parameter is not NULL. The function extracts the raster values in the \code{raster_attribute} parameter for each node, then applies the selected function to this parameter to obtain a single value for each node. For example, mean, sum, modal, min or max. Two of the most popular functions are mean and sum, \code{Default = mean}.
#' @param weighted \code{logical}. If the parameters \bold{raster_attribute and weighted} are \code{TRUE} then the value of the raster attribute for each node is multiplied by its area to obtain an attribute similar to a weighted habitat index.
#' @param area_unit \code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.
#' @param distance A \code{list} of parameters to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. \cr
#'  This list must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. For example: \code{distance(type = "least-cost", resistance = raster_resistance)}. \cr
#' To see more arguments see the \link[Makurhini]{distancefile} function.
#' @param metric A \code{character} indicating the connectivity metric to use: \code{"PC"} (the default and recommended) to calculate the probability of connectivity index, and \code{"IIC"} to calculate the binary integral index of connectivity.
#' @param probability A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the "PC" metric. If \code{probability = NULL}, then a probability of 0.5 will be used.
#' @param distance_thresholds A \code{numeric} indicating the dispersal distance or distances (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range.
#' @param search_buffer \code{numeric}. Distance or distances (i.e., it can be a search distance for each dispersion distance of the parameter \code{distance_thresholds}) used to create a buffer around the focal node (also called focal habitat patch), which is used to select neighbouring nodes (transboundary habitat patches) with which it has the highest probability of connectivity.
#' @param simplify_shape \code{numeric}. It helps to simplify the shape of the focal node by eliminating vertices to buffer and select neighbouring nodes. Its use is recommended when some of the nodes have very complex shapes. See \link[sf]{st_simplify} for details.
#' @param fragmentation \code{logic}. Estimates fragmentation statistics for the focal nodes using the function \link[Makurhini]{MK_Fragmentation}. It is necessary to use the parameters \bold{edge_distance} and \bold{min_node_area}.
#' @param edge_distance  \code{numeric}. Distance to edge in meters. Default equal 500 m. See \link[Makurhini]{MK_Fragmentation} for details.
#' @param min_node_area \code{numeric}. Minimum node area used to calculate the number of nodes with an area smaller than the one provided. Default equal 100 km\out{<sup>2</sup>}. It uses the area units set in the \bold{area_unit} parameter. See \link[Makurhini]{MK_Fragmentation} for details.
#' @param parallel  (\emph{optional, default =} \code{NULL}).
#' A \code{numeric} specifying the number of cores to parallelize the index estimation of the PC or IIC index and its deltas.Particularly useful when you have more than 1000 nodes. By default the analyses are not parallelized.
#' @param save_subfiles \code{character, logical}. Save the result for each focal node in a local folder in \code{.rds} format. If the value is \code{TRUE}, a folder will be generated  in the path specified in the \bold{write} parameter. Otherwise, the folder path must be provided (e.g., \code{'C:/example'}). This parameter is particularly useful in the event of a topological error in the estimation process. For instance, if the function fails due to a topological error in node 601 of 1,000 nodes, the results of the preceding 600 nodes will be saved. You can then correct node 601 and resume your analysis by specifying the address of the folder where your 600 node files were saved.
#' @param write \code{Character} indicating the path and initial prefix of the objects to save, for example, \code{"C:/example/test_focal_"}. By default, nothing is saved. The saved object is a geopackage.
#' @param intern \code{logical}. Show the progress of the process, \code{default = TRUE}. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @references Latorre-Cárdenas, M. C., González-Rodríguez, A., Godínez-Gómez, O., Arima, E. Y., Young, K. R., Denvir, A., ... & Ghilardi, A. (2023). Estimating fragmentation and connectivity patterns of the temperate forest in an avocado-dominated landscape to propose conservation strategies. Land, 12(3), 631.
#' @details
#' A loop is executed, whereby each of the nodes is selected. In each iteration, the following occurs:\cr
#' 1-	When node \code{i} is selected, it becomes a \bold{focal node} (\code{f}).\cr
#' 2-	Then, a buffer is generated with the distance specified in parameter \code{search_buffer}, for example the for example twice the dispersion distance specified in the parameter \code{ distance_thresholds}. This buffer is called \bold{local landscape} and is used to identify neighboring nodes, called \bold{transboundary nodes} (\code{th}).\cr
#' 3-	Next, the index IIC or PC is estimated according to the selected metric using the focal node and the transboundary nodes. This result is referred to as \eqn{\mathbf{IIC}_f} or \eqn{\mathbf{PC}_f}. The index value ranges from 0 to 1, with 1 representing the highest connectivity in the local landscape for the focal node.\cr
#' 4-	Subsequently, the delta \bold{dIIC or dPC} is estimated for the focal node, along with its \code{intra, flux, and connector} deltas.\cr
#' 5-	The function calculates the \bold{Composite Connectivity Index} (\eqn{\mathbf{CCI}_f}) as a prioritization tool for focal nodes. This is based on their individual contribution, weighted by the connectivity in the local landscape:
#' \eqn{\mathbf{CCI}_f = \mathbf{IIC}_f \cdot \mathbf{dIIC}_f} or \eqn{\mathbf{CCI}_f = \mathbf{PC}_f \cdot \mathbf{dPC}_f}. Nodes,  with higher \eqn{\mathbf{CCI}_f} values are found in well-connected local landscapes, making them valuable contributors to connectivity in their immediate landscapes. This makes them ideal candidates for conservation efforts. Conversely, lower \eqn{\mathbf{CCI}_f} values may indicate the need for restoration and conservation actions.
#' 6- In the final step, if the parameter \code{fragmentation} is set to \code{TRUE}, fragmentation statistics are estimated for the local landscape.\cr\cr
#' This process is repeated for each node and stored in an object class sf. For further details, please see Latorre-Cárdenas et al., 2023
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' data("habitat_nodes", package = "Makurhini")
#' nrow(habitat_nodes) # Number of patches
#' test <- MK_Focal_nodes(nodes = habitat_nodes,
#'                          id = "Id",
#'                          attribute = NULL,
#'                          raster_attribute = NULL,
#'                          fun_attribute = NULL,
#'                          distance = list(type = "centroid"),
#'                          metric = "PC",
#'                          probability = 0.5,
#'                          parallel = 4,
#'                          distance_thresholds = 10000,
#'                          search_buffer = 20000,
#'                          fragmentation = TRUE)
#' plot(test["dPC"], breaks = "jenks")
#' plot(test["IComp"], breaks = "jenks")
#' }
#' @note Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv object.size
#' @importFrom purrr map_dbl
#' @importFrom raster raster crop values mask res as.matrix extent raster stack extent<- writeRaster reclassify crs crs<- unique
#' @importFrom terra rast crop values mask vect res
#' @importFrom future multicore multisession plan availableCores
#' @importFrom furrr future_map
#' @importFrom igraph delete_vertices distances
#' @importFrom sf st_bbox st_area st_simplify st_buffer st_intersects st_as_sfc st_as_sf st_geometry st_point_on_surface st_drop_geometry write_sf
#' @importFrom stats median

MK_Focal_nodes <- function(nodes = NULL,
                            id = NULL,
                            attribute = NULL,
                            raster_attribute = NULL,
                            fun_attribute = NULL,
                            weighted = FALSE,
                            area_unit = "ha",
                            distance = list(type= "centroid", resistance = NULL),
                            metric = "IIC",
                            probability = NULL,
                            distance_thresholds = NULL,
                            search_buffer = NULL,
                            simplify_shape = NULL,
                            fragmentation = FALSE,
                            edge_distance = 500,
                            min_node_area = 100,
                            parallel = NULL,
                            write = NULL,
                            save_subfiles = FALSE,
                            intern = TRUE
){
  . = NULL
  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }
  #
  if(is.null(write) & isTRUE(save_subfiles)){
    stop("You need an output path")
  }

  if(is.character(save_subfiles)){
    if(!dir.exists(dirname(save_subfiles))){
      stop("Review save_subfiles path")
    }
  }

  if(is.null(simplify_shape)){
    simplify_shape = 0
  }

  if(any(!(id %in% names(nodes)))){
    stop("Review id parameter")
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

  #Intraconnectivity indicator
  message("1. Nodes attribute")
  if(is.null(attribute)){
    if(!is.null(raster_attribute)){
      if(is.null(fun_attribute)){fun_attribute <- mean}
      if(!is.null(parallel)){
        works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}

        if(.Platform$OS.type == "unix") {
          strat <- future::multicore
        } else {
          strat <- future::multisession
        }

        if(class(raster_attribute)[[1]] == "SpatRaster"){
          raster_attribute <- raster(raster_attribute)
        }

        m <- as.numeric(object.size(raster_attribute))* 0.001

        if(m > 600){
          m <- (m + round(m/3)) *1024^2
          options(future.globals.maxSize = m)
        }

        plan(strategy = strat, gc = TRUE, workers = works)
        indice_tmp <- future_map_dbl(1:nrow(nodes), function(x){
          indice <- raster::crop(raster_attribute, nodes[x,]); indicev <- raster::values(indice)

          if(length(unique(indice)) == 1){
            if(is.na(unique(indice)) | unique(indice) == 0){
              indice3 <- 0
            } else {
              indice3 <- fun_attribute(indicev, na.rm = TRUE)
              if(isTRUE(weighted)){
                indice3 <- indice3 * (unit_convert(st_area(nodes[x,]), "m2", area_unit))
              }
            }
          } else{
            indice <- raster::mask(indice, nodes[x,]); indicev <- raster::values(indice)

            if(fun_attribute == "sum"){
              indice3 <- indicev[!is.na(indicev)] %>%  table(.) %>% as.data.frame(.)
              if(isTRUE(weighted)){
                indice3$resolucion <- unit_convert(raster::res(raster_attribute)[1]^2, "m2", area_unit)
                indice3 <- sum((indice3[[2]]*indice3$resolucion)* as.numeric(as.character(indice3[[1]])))
              } else {
                indice3 <- sum(as.numeric(as.character(indice3[[1]]))*indice3[[2]])
              }
            } else{
              indice3 <- fun_attribute(indicev, na.rm = TRUE)
              if(weighted){
                indice3 <- indice3 * (unit_convert(st_area(nodes[x,]), "m2", area_unit))
              }
            }
          }
          return(indice3)
        }, .progress = intern)

      } else {
        if(class(raster_attribute)[[1]] != "SpatRaster"){
          raster_attribute <- rast(raster_attribute)
        }
        indice_tmp <- map_dbl(1:nrow(nodes), function(x){
          indice <- raster::crop(raster_attribute, nodes[x,]); indicev <- raster::values(indice)

          if(length(unique(indice)) == 1){
            if(is.na(unique(indice)) | unique(indice) == 0){
              indice3 <- 0
            } else {
              indice3 <- fun_attribute(indicev, na.rm = TRUE)
              if(isTRUE(weighted)){
                indice3 <- indice3 * (unit_convert(st_area(nodes[x,]), "m2", area_unit))
              }
            }
          } else{
            indice <- raster::mask(indice, nodes[x,]); indicev <- raster::values(indice)

            if(fun_attribute == "sum"){
              indice3 <- indicev[!is.na(indicev)] %>%  table(.) %>% as.data.frame(.)
              if(isTRUE(weighted)){
                indice3$resolucion <- unit_convert(raster::res(raster_attribute)[1]^2, "m2", area_unit)
                indice3 <- sum((indice3[[2]]*indice3$resolucion)* as.numeric(as.character(indice3[[1]])))
              } else {
                indice3 <- sum(as.numeric(as.character(indice3[[1]]))*indice3[[2]])
              }
            } else{
              indice3 <- fun_attribute(indicev, na.rm = TRUE)
              if(weighted){
                indice3 <- indice3 * (unit_convert(st_area(nodes[x,]), "m2", area_unit))
              }
            }
          }
          return(indice3)
        }, .progress = intern)
      }
      invisible(gc())
    } else {
      indice_tmp <- unit_convert(as.numeric(st_area(nodes)), "m2", "ha")
    }
  } else {
    indice_tmp <- nodes[[attribute]]
  }
  nodes$indice_val <- indice_tmp

  #Conectividad
  message("2. Connectivity indexes")
  if(any(isTRUE(save_subfiles)|is.character(save_subfiles))){
    '%!in%' <- function(x,y)!('%in%'(x,y))
    if(is.character(save_subfiles)){
      if(!dir.exists(save_subfiles)){
        dir.create(save_subfiles)
      }
      ss <- save_subfiles; iter <- list.files(ss, pattern = ".Rds$")
      if(length(iter) > 0){
        salida <- if(is.character(write)){basename(write)} else {basename(ss)}
        iter <- gsub(salida, "", iter) %>% gsub("_|.Rds|MK", "", .) %>% as.numeric(.)
        iter2 <- which(1:nrow(nodes) %!in% iter); iter <- 1:nrow(nodes); iter <- iter[iter2]
        if(length(iter) == 0){
          iter <- nrow(nodes)
        }
      } else {
        iter <- 1:nrow(nodes)
      }

    } else {
      dname <- dir(dirname(write), pattern = basename(write))
      if(length(dname) > 0){
        dname <- gsub(paste0(basename(write),"|_"), "", dname) %>% as.numeric() %>% max(.)
        dname <- paste0("_", dname + 1) %>% gsub("__", "_", .); ss <- paste0(write, dname); dir.create(ss)
      } else {
        ss <- dirname(write)
        dname <- paste0(basename(write),"_", 1); ss <- paste0(ss, "/", gsub("__", "_", dname)); dir.create(ss)
      }
      iter <- 1:nrow(nodes)
    }
  } else {
    iter <- 1:nrow(nodes)
  }

  if(length(iter) > nrow(nodes)){
    iter <- 1:nrow(nodes)
  }

  if(is.null(parallel)){
    resultado <- tryCatch(lapply(iter, function(i){
      if(isTRUE(intern)){
        cat(paste0("\r", i, " of ", nrow(nodes), " nodes processed"))
      }
      i.0 <- TopoClean(nodes[i,], xsimplify = TRUE)
      AI <- st_buffer(i.0, max(search_buffer)); i.1 <- nodes[as.vector(st_intersects(AI, nodes, sparse = FALSE)),]
      if(simplify_shape > 0){
        i.1 <- st_simplify(i.1, dTolerance = simplify_shape, preserveTopology = TRUE)
      }
      i.1$gridcode <- NULL

      if(!is.null(distance$resistance)){
        distance1 <- distance; distance1$resistance <- crop(distance$resistance, st_bbox(i.1) %>% st_as_sfc() %>%
                                                              st_as_sf() %>% st_buffer(max(distance_thresholds)))
      } else {
        distance1 <- distance; distance1$resistance <- NULL
      }
      if(nrow(i.1) > 1){
        puntos <- data.frame(Id = i.1[[id]], lon = map_dbl(st_geometry(i.1), ~st_point_on_surface(.x)[[1]]),
                             lat = map_dbl(st_geometry(i.1), ~st_point_on_surface(.x)[[2]])) %>%
          st_as_sf(coords = c("lon", "lat"), crs = st_crs(i.1), agr = "constant")
        distancias <- distancefile(nodes = puntos,
                                   id = "Id",
                                   type = distance1$type,
                                   keep = distance$keep,
                                   resistance = distance1$resistance,
                                   resist.units = distance$resist.units,
                                   CostFun = distance$CostFun,
                                   ngh = distance$ngh,
                                   mask = distance$mask,
                                   threshold = distance$threshold,
                                   geometry_out = distance$geometry_out,
                                   bounding_circles = distance$bounding_circles,
                                   parallel = distance$parallel,
                                   ActiveParallel = FALSE,
                                   least_cost.java = distance1$least_cost.java,
                                   cores.java = distance1$cores.java,
                                   ram.java = distance1$ram.java,
                                   pairwise = FALSE,
                                   write = NULL)
      }
      AI <- i.0; rm(i.0)
      resultado.i <- lapply(1:length(distance_thresholds), function(j){
        if(search_buffer[j] == max(search_buffer)){
          j.1 <- i.1; j.1$gridcode <- NULL
        } else {
          j.1 <- i.1[as.vector(st_intersects(st_buffer(AI, search_buffer[j]), i.1, sparse = FALSE)),]; j.1$gridcode <- NULL
        }

        if(nrow(j.1) > 2){
          if(is.null(raster_attribute)){
            LA <- st_bbox(j.1) %>% st_as_sfc() %>% st_as_sf() %>% st_area() %>% as.numeric(); LA <- unit_convert(LA, "m2", "ha")
          } else {
            LA <- st_bbox(j.1) %>% st_as_sfc() %>% st_as_sf() %>% st_buffer(res(raster_attribute)[1]); LA <- mask(raster_attribute, LA)
            LA <- freq(LA) %>% as.data.frame(); LA <- sum(LA$count[which(!is.na(LA$value))])
          }

          distance2 <- distancias[which(row.names(distancias) %in% as.character(j.1[[id]])),which(colnames(distancias) %in% as.character(j.1[[id]]))]
          rownames(distance2) <- NULL; colnames(distance2) <- NULL

          dIIC <- tryCatch(MK_dPCIIC(nodes = j.1,
                                     attribute = "indice_val",
                                     area_unit = "ha",
                                     distance = distance2,
                                     parallel = NULL,
                                     overall = TRUE,
                                     onlyoverall = FALSE,
                                     LA = LA,
                                     metric = metric,
                                     probability = probability,
                                     distance_thresholds = distance_thresholds[j],
                                     id_sel = which(j.1[[id]] == nodes[[id]][i]),
                                     intern = FALSE),
                           error = function(err) err)

          if(inherits(dIIC, "error")){
            if(nrow(j.1)>1000){
              a <- ps::ps()
              a <- a[which(a$name =="Rscript.exe"),1]

              for(q in a){
                tools::pskill(q)
              }
            }

            j.1 <- st_simplify(j.1, dTolerance = simplify_shape, preserveTopology = TRUE)

            if(!is.null(distance$resistance)){
              distance2 <- distance1; distance2$resistance <- crop(distance1$resistance, st_bbox(j.1) %>% st_as_sfc() %>%
                                                                     st_as_sf() %>% st_buffer(distance_thresholds[j]))
            }

            dIIC <- tryCatch(MK_dPCIIC(nodes = j.1,
                                       attribute = "indice_val",
                                       area_unit = "ha",
                                       distance = distance2,
                                       parallel = NULL,
                                       overall = TRUE,
                                       onlyoverall = FALSE,
                                       LA = LA,
                                       metric = metric,
                                       probability = probability,
                                       distance_thresholds = distance_thresholds[j],
                                       id_sel = which(j.1[[id]] == nodes[[id]][i]),
                                       intern = FALSE),
                             error = function(err) err)

            if(inherits(dIIC, "error")){
              if(nrow(j.1)>1000){
                a <- ps::ps()
                a <- a[which(a$name =="Rscript.exe"),1]

                for(q in a){
                  tools::pskill(q)
                }
              }
              stop(paste0("checa parche en la fila ", i))
            } else {
              #IIC
              IIC.2 <- dIIC[[2]] %>% as.data.frame(); IIC <- t(IIC.2[,2]) %>% as.data.frame()
              names(IIC) <- IIC.2[,1]; IIC <- IIC[,2:3]

              #dIIC
              dIIC <- dIIC[[1]] %>% st_drop_geometry()
              dIIC <- dIIC[which(dIIC[[grep(id, names(dIIC))]] == nodes[[id]][j]),]
              dIIC <- dIIC[,c((ncol(dIIC)-3):ncol(dIIC))]; dIIC <- cbind(IIC, dIIC)
            }
          } else {
            #IIC
            IIC.2 <- dIIC[[2]] %>% as.data.frame(); IIC <- t(IIC.2[,2]) %>% as.data.frame()
            names(IIC) <- IIC.2[,1]; IIC <- IIC[,2:3]

            #dIIC
            dIIC <- dIIC[[1]] %>% st_drop_geometry(); dIIC <- dIIC[which(dIIC[[grep(id, names(dIIC))]] == nodes[[id]][i]),]
            dIIC <- dIIC[,c((ncol(dIIC)-3):ncol(dIIC))]; dIIC <- cbind(IIC, dIIC)
          }

          if(isTRUE(fragmentation)){
            x.2 <- MK_Fragmentation(nodes = j.1, edge_distance = edge_distance,
                                    min_node_area = min_node_area, area_unit = area_unit,
                                    perimeter_unit = "km", plot = FALSE, write = NULL)
            x.2 <- t(as.data.frame(x.2[[1]])) |> as.data.frame(x = _)
            names(x.2) <- as.character(x.2[1,]); x.2 <- x.2[2,]
            x.2[] <- lapply(x.2, function(y) {as.numeric(as.character(y))})
            dIIC <- cbind(dIIC, x.2)
          }
        } else {
          IIC <- data.frame('EC(IIC)' = NA, IIC =  NA, check.names = FALSE)
          dIIC <- data.frame(dIIC = NA, dIICintra = NA, dIICflux = NA, dIICconnector = NA,
                             check.names = FALSE)
          if(metric == "PC"){
            names(IIC) <- c('EC(PC)', "PC"); names(dIIC) <- gsub("IIC", "PC", names(dIIC))
          }
          dIIC <- cbind(IIC, dIIC)
          if(isTRUE(fragmentation)){
            if(nrow(j.1) == 1){
              x.2 <- MK_Fragmentation(nodes = j.1, edge_distance = edge_distance,
                                      min_node_area = min_node_area, area_unit = area_unit,
                                      perimeter_unit = "km", plot = FALSE, write = NULL)
              x.2 <- t(as.data.frame(x.2[[1]])) |> as.data.frame(x = _)
              names(x.2) <- as.character(x.2[1,]); x.2 <- x.2[2,]
              x.2[] <- lapply(x.2, function(y) {as.numeric(as.character(y))})
              dIIC <- cbind(dIIC, x.2)
            } else {
              x.2 <- data.frame("Patch area (ha)" = NA, "Number of patches" = NA,
                                "Size (mean)" = NA, "Patches < minimum patch area" = NA,
                                "Patches < minimum patch area (%)" = NA, "Total edge" = NA,
                                "Edge density" = NA, "Patch density" = NA, "Total Core Area (ha)" = NA,
                                "Cority" = NA, "Shape Index (mean)" = NA,
                                "FRAC (mean)" = NA, "MESH (ha)" = NA, check.names = FALSE)
              dIIC <- cbind(dIIC, x.2)
            }
          }
        }
        return(dIIC)
      })

      if(any(isTRUE(save_subfiles) | is.character(save_subfiles))){
        salida <- paste0(if(is.character(save_subfiles)){basename(ss)} else {basename(write)},
                         "_MK_", i,".Rds") |> gsub("__", "_", x = _)
        salida <- paste0(ss, "/", salida); saveRDS(resultado.i, file = salida)
      }
      return(resultado.i)
    }), error = function(err)err)
  } else {
    if(is.null(raster_attribute)){
      works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}

      if(.Platform$OS.type == "unix") {
        strat <- future::multicore
      } else {
        strat <- future::multisession
      }

      if(class(raster_attribute)[[1]] == "SpatRaster"){
        raster_attribute <- raster(raster_attribute)
      }

      m <- as.numeric(object.size(raster_attribute))* 0.001

      if(m > 600){
        m <- (m + round(m/3)) *1024^2
        options(future.globals.maxSize = m)
      }

      plan(strategy = strat, gc = TRUE, workers = works)
    }
    resultado <- tryCatch(future_map(iter, function(i){
      i.0 <- TopoClean(nodes[i,], xsimplify = TRUE)
      AI <- st_buffer(i.0, max(search_buffer)); i.1 <- nodes[as.vector(st_intersects(AI, nodes, sparse = FALSE)),]
      if(simplify_shape > 0){
        i.1 <- st_simplify(i.1, dTolerance = simplify_shape, preserveTopology = TRUE)
      }
      i.1$gridcode <- NULL

      if(!is.null(distance$resistance)){
        distance1 <- distance; distance1$resistance <- crop(distance$resistance, st_bbox(i.1) %>% st_as_sfc() %>%
                                                              st_as_sf() %>% st_buffer(max(distance_thresholds)))
      } else {
        distance1 <- distance; distance1$resistance <- NULL
      }
      if(nrow(i.1) > 1){
        puntos <- data.frame(Id = i.1[[id]], lon = map_dbl(st_geometry(i.1), ~st_point_on_surface(.x)[[1]]),
                             lat = map_dbl(st_geometry(i.1), ~st_point_on_surface(.x)[[2]])) %>%
          st_as_sf(coords = c("lon", "lat"), crs = st_crs(i.1), agr = "constant")
        distancias <- distancefile(nodes = puntos,
                                   id = "Id",
                                   type = distance1$type,
                                   keep = distance$keep,
                                   resistance = distance1$resistance,
                                   resist.units = distance$resist.units,
                                   CostFun = distance$CostFun,
                                   ngh = distance$ngh,
                                   mask = distance$mask,
                                   threshold = distance$threshold,
                                   geometry_out = distance$geometry_out,
                                   bounding_circles = distance$bounding_circles,
                                   parallel = NULL,
                                   ActiveParallel = FALSE,
                                   least_cost.java = distance1$least_cost.java,
                                   cores.java = distance1$cores.java,
                                   ram.java = distance1$ram.java,
                                   pairwise = FALSE,
                                   write = NULL)
      }
      AI <- i.0; rm(i.0)
      resultado.i <- lapply(1:length(distance_thresholds), function(j){
        if(search_buffer[j] == max(search_buffer)){
          j.1 <- i.1; j.1$gridcode <- NULL
        } else {
          j.1 <- i.1[as.vector(st_intersects(st_buffer(AI, search_buffer[j]), i.1, sparse = FALSE)),]; j.1$gridcode <- NULL
        }

        if(nrow(j.1) > 2){
          if(is.null(raster_attribute)){
            LA <- st_bbox(j.1) %>% st_as_sfc() %>% st_as_sf() %>% st_area() %>% as.numeric(); LA <- unit_convert(LA, "m2", "ha")
          } else {
            LA <- st_bbox(j.1) %>% st_as_sfc() %>% st_as_sf() %>% st_buffer(res(raster_attribute)[1]); LA <- mask(raster_attribute, LA)
            LA <- freq(LA) %>% as.data.frame(); LA <- sum(LA$count[which(!is.na(LA$value))])
          }

          distance2 <- distancias[which(row.names(distancias) %in% as.character(j.1[[id]])),which(colnames(distancias) %in% as.character(j.1[[id]]))]
          rownames(distance2) <- NULL; colnames(distance2) <- NULL

          dIIC <- tryCatch(MK_dPCIIC(nodes = j.1,
                                     attribute = "indice_val",
                                     area_unit = "ha",
                                     distance = distance2,
                                     parallel = NULL,
                                     overall = TRUE,
                                     onlyoverall = FALSE,
                                     LA = LA,
                                     metric = metric,
                                     probability = probability,
                                     distance_thresholds = distance_thresholds[j],
                                     id_sel = which(j.1[[id]] == nodes[[id]][i]),
                                     intern = FALSE),
                           error = function(err) err)

          if(inherits(dIIC, "error")){
            if(nrow(j.1)>1000){
              a <- ps::ps()
              a <- a[which(a$name =="Rscript.exe"),1]

              for(q in a){
                tools::pskill(q)
              }
            }

            j.1 <- st_simplify(j.1, dTolerance = simplify_shape, preserveTopology = TRUE)

            if(!is.null(distance$resistance)){
              distance2 <- distance1; distance2$resistance <- crop(distance1$resistance, st_bbox(j.1) %>% st_as_sfc() %>%
                                                                     st_as_sf() %>% st_buffer(distance_thresholds[j]))
            }

            dIIC <- tryCatch(MK_dPCIIC(nodes = j.1,
                                       attribute = "indice_val",
                                       area_unit = "ha",
                                       distance = distance2,
                                       parallel = NULL,
                                       overall = TRUE,
                                       onlyoverall = FALSE,
                                       LA = LA,
                                       metric = metric,
                                       probability = probability,
                                       distance_thresholds = distance_thresholds[j],
                                       id_sel = which(j.1[[id]] == nodes[[id]][i]),
                                       intern = FALSE),
                             error = function(err) err)

            if(inherits(dIIC, "error")){
              if(nrow(j.1)>1000){
                a <- ps::ps()
                a <- a[which(a$name =="Rscript.exe"),1]

                for(q in a){
                  tools::pskill(q)
                }
              }
              stop(paste0("checa parche en la fila ", i))
            } else {
              #IIC
              IIC.2 <- dIIC[[2]] %>% as.data.frame(); IIC <- t(IIC.2[,2]) %>% as.data.frame()
              names(IIC) <- IIC.2[,1]; IIC <- IIC[,2:3]

              #dIIC
              dIIC <- dIIC[[1]] %>% st_drop_geometry()
              dIIC <- dIIC[which(dIIC[[grep(id, names(dIIC))]] == nodes[[id]][j]),]
              dIIC <- dIIC[,c((ncol(dIIC)-3):ncol(dIIC))]; dIIC <- cbind(IIC, dIIC)
            }
          } else {
            #IIC
            IIC.2 <- dIIC[[2]] %>% as.data.frame(); IIC <- t(IIC.2[,2]) %>% as.data.frame()
            names(IIC) <- IIC.2[,1]; IIC <- IIC[,2:3]

            #dIIC
            dIIC <- dIIC[[1]] %>% st_drop_geometry(); dIIC <- dIIC[which(dIIC[[grep(id, names(dIIC))]] == nodes[[id]][i]),]
            dIIC <- dIIC[,c((ncol(dIIC)-3):ncol(dIIC))]; dIIC <- cbind(IIC, dIIC)
          }

          if(isTRUE(fragmentation)){
            x.2 <- MK_Fragmentation(nodes = j.1, edge_distance = edge_distance,
                                    min_node_area = min_node_area, area_unit = area_unit,
                                    perimeter_unit = "km", plot = FALSE, write = NULL)
            x.2 <- t(as.data.frame(x.2[[1]])) |> as.data.frame(x = _)
            names(x.2) <- as.character(x.2[1,]); x.2 <- x.2[2,]
            x.2[] <- lapply(x.2, function(y) {as.numeric(as.character(y))})
            dIIC <- cbind(dIIC, x.2)
          }
        } else {
          IIC <- data.frame('EC(IIC)' = NA, IIC =  NA, check.names = FALSE)
          dIIC <- data.frame(dIIC = NA, dIICintra = NA, dIICflux = NA, dIICconnector = NA,
                             check.names = FALSE)
          if(metric == "PC"){
            names(IIC) <- c('EC(PC)', "PC"); names(dIIC) <- gsub("IIC", "PC", names(dIIC))
          }
          dIIC <- cbind(IIC, dIIC)
          if(isTRUE(fragmentation)){
            if(nrow(j.1) == 1){
              x.2 <- MK_Fragmentation(nodes = j.1, edge_distance = edge_distance,
                                      min_node_area = min_node_area, area_unit = area_unit,
                                      perimeter_unit = "km", plot = FALSE, write = NULL)
              x.2 <- t(as.data.frame(x.2[[1]])) |> as.data.frame(x = _)
              names(x.2) <- as.character(x.2[1,]); x.2 <- x.2[2,]
              x.2[] <- lapply(x.2, function(y) {as.numeric(as.character(y))})
              dIIC <- cbind(dIIC, x.2)
            } else {
              x.2 <- data.frame("Patch area (ha)" = NA, "Number of patches" = NA,
                                "Size (mean)" = NA, "Patches < minimum patch area" = NA,
                                "Patches < minimum patch area (%)" = NA, "Total edge" = NA,
                                "Edge density" = NA, "Patch density" = NA, "Total Core Area (ha)" = NA,
                                "Cority" = NA, "Shape Index (mean)" = NA,
                                "FRAC (mean)" = NA, "MESH (ha)" = NA, check.names = FALSE)
              dIIC <- cbind(dIIC, x.2)
            }
          }
        }
        return(dIIC)
      })

      if(any(isTRUE(save_subfiles) | is.character(save_subfiles))){
        salida <- paste0(if(is.character(save_subfiles)){basename(ss)} else {basename(write)},
                         "_MK_", i,".Rds") |> gsub("__", "_", x = _)
        salida <- paste0(ss, "/", salida); saveRDS(resultado.i, file = salida)
      }
      return(resultado.i)
    }, .progress = intern), error = function(err)err)
    close_multiprocess(works)
  }

  if(inherits(resultado, "error")){
    stop(resultado)
  } else {
    if(any(isTRUE(save_subfiles)|is.character(save_subfiles))){
      salida <- paste0(if(is.character(save_subfiles)){basename(ss)} else {basename(write)})
      resu <- list.files(ss, pattern = ".Rds$", full.names = TRUE)
      resu.2 <- basename(resu); resu.2 <- gsub(paste0(".Rds|", salida, "_MK_"), "", resu.2)
      resultado <-  lapply(1:length(resu), function(x){
        x.0 <- resu[which(resu.2 == as.character(x))]
        x.1 <- readRDS(x.0); return(x.1)
      })
    }

    if(length(distance_thresholds) > 1){
      IIC_dIIC <- lapply(1:length(distance_thresholds), function(x){
        x.1 <- cbind(nodes, map_dfr(resultado, function(y){return(y[[x]])}))
        rownames(x.1) <- NULL
        x.1$IComp <- x.1[[grep("EC", names(x.1))+1]] * x.1[[grep("intra", names(x.1))-1]]
        select_column <- if(isTRUE(fragmentation)){
          c(1:grep("connector", names(x.1)),
            ncol(x.1), (grep("connector", names(x.1))+1):(ncol(x.1)-2))
        } else {
          c(1:(ncol(x.1)-2), ncol(x.1))
        }
        x.1 <- x.1[,select_column]
        if(!is.null(write)){
          salida <- paste0(write, "_d", distance_thresholds[x], ".gpkg") |> gsub("__", "_", x = _)
          write_sf(IIC_dIIC, salida)
        }
        return(x.1)
      })
      message("Connectivity indexes ready!")
      names(IIC_dIIC) <- paste0("d", distance_thresholds)
      if(!is.null(write)){
        salida <- paste0(write, "_", paste0(distance_thresholds, collapse = "_"), ".rds") |> gsub("__", "_", x = _)
        saveRDS(IIC_dIIC, file = salida)
      }
    } else {
      IIC_dIIC <- cbind(nodes, map_dfr(resultado, function(x){x[[1]]}))
      rownames(IIC_dIIC) <- NULL; names(IIC_dIIC) <- gsub("\\.\\.|\\.\\.\\.|\\.\\.\\.\\.", ".", names(IIC_dIIC))
      #Indice compuesto
      IIC_dIIC$IComp <- IIC_dIIC[[grep("EC", names(IIC_dIIC))+1]] * IIC_dIIC[[grep("intra", names(IIC_dIIC))-1]]
      message("Indices de conectividad listos")
      select_column <- if(isTRUE(fragmentation)){
        c(1:grep("connector", names(IIC_dIIC)),
          ncol(IIC_dIIC),
          (grep("connector", names(IIC_dIIC))+1):(ncol(IIC_dIIC)-2))
      } else {
        c(1:(ncol(IIC_dIIC)-2),ncol(IIC_dIIC))
      }
      IIC_dIIC <- IIC_dIIC[,select_column]
      if(!is.null(write)){
        salida <- paste0(write, "_", distance_thresholds, ".gpkg") |> gsub("__", "_", x = _)
        write_sf(IIC_dIIC, salida)
      }
    }
  }
  if(isTRUE(intern)){
    message("Done!")
  }
  return(IIC_dIIC)
}

