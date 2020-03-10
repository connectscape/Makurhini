#' Estimate the integral index of connectivity (IIC) or the probability of connectivity (PC)
#'
#' Use this function to calculate the PC and IIC indexes under one or several distance thresholds.
#' @param nodes Object of class sf, SpatialPolygonsDataFrame or raster.
#' It must be in a projected coordinate system.
#' @param id character. If nodes is a shapefile then you can specify the column name with the node ID.
#'  If NULL, then a new id will be generated. If nodes is a raster layer then raster values (Integer) will be taken as "id".
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute in the data table selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit (e.g., "km2", "cm2", "ha";
#'  see Makurhini::unit_convert). Default equal to hectares "ha".
#' @param restauration character or vector. If nodes is a shappefile then you must specify the name of the column
#' with restauration value. If nodes is a raster layer then must be a numeric vector with restauration values
#' to each node in the raster. Binary values (0,1), where 1 = existing nodes in the landscape, and 0 = a new node
#'  to add to the initial landscape (restored).
#' @param distance list. Distance parameters. For example: type, resistance,or keep. For "type" choose one
#' of the distances: "centroid" (faster), "edge", "hausdorff-edge",
#' "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have
#' to use the "resistance" argument. To see more arguments see the ?distancefile.
#' @param metric character. Choose a connectivity metric: "IIC" considering topologycal distances or "PC" considering maximum product probabilities.
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param overall logical. If TRUE, then the EC index will be added to the result which is transformed into a list. Default equal to FALSE
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to ha).
#' @param dA logical. If TRUE, then the delta attribute will be added to the node's importance result.
#' @param dvars logical. If TRUE, then the absolute variation will be added to the node's importance result.
#' @param parallel logical. Parallelize the function using furrr package and multiprocess plan, default = FALSE.
#' @param rasterparallel logical. If parallel is FALSE and nodes is a raster then you can use this argument to assign the metrics values to the nodes raster. It is useful when raster resolution is less than 100 m2.
#' @param write character. Write output shapefile and overall table (if TRUE overall argument).
#'   It is necessary to specify the "Folder direction" + "Initial prefix",  for example, "C:/ejemplo".
#' @references Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid. Available at \url{www.conefor.org}.\cr
#' Pascual-Hortal, L. & Saura, S. 2006. Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology 21 (7): 959-967.\cr
#' Saura, S. & Pascual-Hortal, L. 2007. A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning 83 (2-3): 91-103.
#' @export
#' @examples
#' library(Makurhini)
#' data("vegetation_patches", package = "Makurhini")
#' nrow(vegetation_patches) # Number of patches
#' #Two distance threshold,
#' IIC <- MK_dPCIIC(nodes = vegetation_patches, id = "id", attribute = NULL,
#'                     distance = list(type = "centroid"),
#'                     metric = "IIC", distance_thresholds = c(5000, 10000)) # 5 and 10 km
#' IIC
#' #Using raster
#' data("raster_vegetation_patches", package = "Makurhini")
#' PC <- MK_dPCIIC(nodes = raster_vegetation_patches, id = NULL, attribute = NULL,
#'                     distance = list(type = "centroid"),
#'                     metric = "PC", probability = 0.5,
#'                     distance_thresholds = 40000) # 40 km
#' PC
#'
#' PC_parallel <- MK_dPCIIC(nodes = raster_vegetation_patches, id = NULL, attribute = NULL,
#'                     distance = list(type = "centroid"),
#'                     metric = "PC", probability = 0.5,
#'                     distance_thresholds = c(40000, 60000), parallel=TRUE) #40 and 60 km
#' PC_parallel
#' @importFrom dplyr progress_estimated
#' @importFrom methods as
#' @importFrom utils write.table warnErrList
#' @importFrom iterators iter
#' @importFrom foreach foreach %dopar%
#' @importFrom purrr map
#' @importFrom raster values as.matrix extent raster stack extent<- writeRaster reclassify crs crs<-
#' @importFrom future multiprocess plan
#' @importFrom furrr future_map
MK_dPCIIC <- function(nodes, id = NULL, attribute  = NULL,
                      area_unit = "ha", restauration = NULL,
                      distance = list(type= "centroid", resistance = NULL),
                      metric = c("IIC", "PC"),
                      probability, distance_thresholds = NULL,
                      overall = FALSE, dA = FALSE, dvars =FALSE,
                      LA = NULL, parallel = FALSE, rasterparallel = FALSE, write = NULL) {
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
    if (is.null(probability) | !is.numeric(probability)) {
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

  if(class(nodes)[1] == "sf") {
    nodes <- as(nodes, 'Spatial')
  }

  options(warn = -1)
  ttt.2 <- getwd()
  #

  temp.1 <- paste0(tempdir(), "/TempInputs", sample(1:1000, 1, replace = TRUE))
  dir.create(temp.1, recursive = T)

  #
  if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
    if (is.null(id)) {
      nodes@data$IdTemp <- 1:nrow(nodes)
    } else {
      nodes@data$IdTemp <- nodes@data[[id]]
    }
    id = "IdTemp"
  } else {
    id = NULL
  }

  nodesfile(nodes, id = id, attribute = attribute, area_unit = area_unit,
            multiple = NULL, restauration = restauration,
            prefix=NULL, write = paste0(temp.1,"/nodes.txt"))

  distancefile(nodes,  id = id, type = distance$type, keep = distance$keep,
               resistance = distance$resistance, CostFun = distance$CostFun, ngh = distance$ngh,
               threshold = distance$threshold, mask = distance$mask,
               distance_unit = distance$distance_unit, geometry_out = distance$geometry_out,
               write = paste0(temp.1,"/Dist.txt"))

  setwd(temp.1)
  if(!is.null(restauration)){
    rest = TRUE
  } else {
    rest = FALSE
  }
  #
  if (is.null(distance$threshold)) {
    pairs = "all"
  } else {
    pairs = "notall"
  }

  x = NULL
  if(isFALSE(parallel)){
    pb <- progress_estimated(length(distance_thresholds), 0)
    result <- foreach(x=iter(distance_thresholds), .errorhandling = 'pass') %dopar%
      {
        if(length(distance_thresholds)>1){
          pb$tick()$print()
        }

        dMetric <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                              typeconnection = "dist", typepairs = pairs,
                              index = metric, thdist = x,
                              multdist = NULL, conprob = probability,
                              onlyoverall = FALSE, LA = LA,
                              nrestauration = rest,
                              prefix = NULL, write = NULL)

        if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
          Merge_metric <- merge_conefor(datat = dMetric[[which(map(dMetric, function(x) ncol(x)) >= 11)]], pattern = NULL,
                                        merge_shape = nodes, id = "IdTemp",
                                        write = if (!is.null(write)) paste0(write, "_d", x,".shp"),
                                        dA = dA, var = dvars)

          Merge_metric$"IdTemp" <- NULL

        } else {
          datat <- dMetric[[which(map(dMetric, function(x) ncol(x)) >= 11)]] %>% as.data.frame(.)
          datat <- datat[, as.numeric(which(colSums(!is.na(datat)) > 0))]

          if(isFALSE(dA)){
            datat$dA <- NULL
            datat$varA <- NULL
            datat[,which(unique(is.na(datat)) == TRUE)] <- NULL
          }

          if (isFALSE(dvars)){
            datat <- datat[which(grepl("var",names(datat), fixed = TRUE) == FALSE)]
          }

          ###rasters values
          rp <- unique(raster::values(nodes))
          rp <- as.vector(rp)
          rp <- rp[which(!is.na(rp))]

          if(isTRUE(rasterparallel)){
            m <- matrix(nrow = nrow(datat), ncol = 2)
            m[,1] <- datat[,1]
            plan(strategy = multiprocess)
            r_metric <- future_map(2:ncol(datat), function(c){
              x1 <- datat[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)}, .progress = TRUE)
            future:::ClusterRegistry("stop")
          } else {
            m <- matrix(nrow = nrow(datat), ncol = 2)
            m[,1] <- datat[,1]
            r_metric <- map(2:ncol(datat), function(c){
              x1 <- datat[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)})
          }

          r_metric_2 <- list()

          for(i in 2:(length(r_metric)+1)){
            r_metric_2[[i]] <- r_metric[[i-1]]
          }

          r_metric_2[[1]] <- nodes
          r_metric_2 <- stack(r_metric_2)
          names(r_metric_2) <- names(datat[,1:ncol(datat)])

          if (!is.null(write)){
            n <- names(datat[,1:ncol(datat)])
            n <- map(as.list(2:length(n)), function(w){
              x1 <- r_metric_2[[w]]
              crs(x1) <- crs(r_metric_2)
              writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x, ".tif"), overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
            })
          }

        }

        if(isTRUE(overall)){
          roverall <- dMetric[[which(map(dMetric, function(x) paste0(nrow(x), ncol(x))) == "32" | map(dMetric, function(x) paste0(nrow(x), ncol(x))) == "22")]]
          names(roverall) <- c("Index", "Value")

          result_interm <- list()

          if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
            result_interm[[1]] <- Merge_metric
          } else {
            result_interm[[1]] <- r_metric_2
          }
          result_interm[[2]] <- roverall
          names(result_interm) <- c(paste0("_node_importances_d",x), paste0("_overall_d", x))

          if (!is.null(write)){
            write.table(roverall, paste0(write, "_overall_d", x), sep="\t",
                        row.names = FALSE, col.names = TRUE)
          }
        } else {
          if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
            result_interm <- Merge_metric
          } else {
            result_interm <- r_metric_2
          }
        }
        return(result_interm)
      }

  } else {
    plan(strategy = multiprocess, gc = TRUE)
    result <- future_map(distance_thresholds, function(x){

      dMetric <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                            typeconnection = "dist", typepairs = pairs,
                            index = metric, thdist = x,
                            multdist = NULL, conprob = probability,
                            onlyoverall = FALSE, LA = LA,
                            nrestauration = rest,
                            prefix = NULL, write = NULL)

      if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
        Merge_metric <- merge_conefor(datat = dMetric[[which(map(dMetric, function(x) ncol(x)) >= 11)]], pattern = NULL,
                                      merge_shape = nodes, id = "IdTemp",
                                      write = if (!is.null(write)) paste0(write, "_d", x,".shp"),
                                      dA = dA, var = dvars)

        Merge_metric$"IdTemp" <- NULL

      } else {
        datat <- dMetric[[which(map(dMetric, function(x) ncol(x)) >= 11)]] %>% as.data.frame(.)
        datat <- datat[, as.numeric(which(colSums(!is.na(datat)) > 0))]

        if(isFALSE(dA)){
          datat$dA <- NULL
          datat$varA <- NULL
          datat[,which(unique(is.na(datat)) == TRUE)] <- NULL
        }

        if (isFALSE(dvars)){
          datat <- datat[which(grepl("var",names(datat), fixed = TRUE) == FALSE)]
        }

        ###rasters values
        rp <- unique(raster::values(nodes))
        rp <- as.vector(rp)
        rp <- rp[which(!is.na(rp))]

        m <- matrix(nrow = nrow(datat), ncol = 2)
        m[,1] <- datat[,1]
        r_metric <- map(2:ncol(datat), function(c){
          x1 <- datat[,c(1, c)]
          for(i in rp){
            m[which(m == i),2] <- x1[which(x1[,1]== i),2]
            }
          x1 <- reclassify(nodes, rcl = m)
          return(x1)})

        r_metric_2 <- list()

        for(i in 2:(length(r_metric)+1)){
          r_metric_2[[i]] <- r_metric[[i-1]]
        }

        r_metric_2[[1]] <- nodes
        r_metric_2 <- stack(r_metric_2)
        names(r_metric_2) <- names(datat[,1:ncol(datat)])

        if (!is.null(write)){
          n <- names(datat[,1:ncol(datat)])
          n <- map(as.list(2:length(n)), function(w){
            x1 <- r_metric_2[[w]]
            crs(x1) <- crs(r_metric_2)
            writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x, ".tif"), overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
          })
        }

      }

      if(isTRUE(overall)){
        roverall <- dMetric[[which(map(dMetric, function(x) paste0(nrow(x), ncol(x))) == "32" | map(dMetric, function(x) paste0(nrow(x), ncol(x))) == "22")]]
        names(roverall) <- c("Index", "Value")

        result_interm <- list()

        if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
          result_interm[[1]] <- Merge_metric
        } else {
          result_interm[[1]] <- r_metric_2
        }
        result_interm[[2]] <- roverall
        names(result_interm) <- c(paste0("_node_importances_d",x), paste0("_overall_d", x))

        if (!is.null(write)){
          write.table(roverall, paste0(write, "_overall_d", x), sep="\t",
                      row.names = FALSE, col.names = TRUE)
        }
      } else {
        if(class(nodes)[1] == "SpatialPolygonsDataFrame"){
          result_interm <- Merge_metric
        } else {
          result_interm <- r_metric_2
        }
      }
      return(result_interm)
    }, .progress = TRUE)
    future:::ClusterRegistry("stop")
  }

  #
  if(!is.null(attr(warnErrList(result), "warningMsg")[[1]])) {
    setwd(ttt.2)
    stop(warnErrList(result))
  } else {
    if (isTRUE(isFALSE(overall) && length(distance_thresholds) == 1)){
      result <- result[[1]]
    } else {
      names(result) <- paste0("d", distance_thresholds)
    }
    setwd(ttt.2)
  }
  return(result)
}
