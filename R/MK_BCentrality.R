#' Betweenness Centrality metrics using Conefor command line
#'
#' Use this function to calculate the BC, BCIIC and BCPC indexes under one or several distance thresholds.
#' @param nodes Object of class sf, SpatialPolygonsDataFrame or raster.
#'  The shapefile must be in a projected coordinate system.
#' @param id character. If nodes is a shapefile then you can specify the column name with the node ID.
#'  If NULL, then a new id will be generated.
#' @param attribute character or vector. If nodes is a shappefile then you must specify the column name
#'  with the attribute in the data table selected for the nodes. If nodes is a raster layer then it must be
#'  a numeric vector with the node's attribute. The length of the vector must be equal to the number of nodes.
#'   The numeric vector is multiplied by the area of each node to obtain a weighted habitat index.
#'   If NULL the node area will be used as a node attribute, the unit area can be selected using the "area_unit" argument.
#' @param area_unit character. If attribute is NULL you can set an area unit, ?unit_covert
#' compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to "m2".
#' @param distance list. Distance parameters. For example: type, resistance,or keep. For "type" choose one
#' of the distances: "centroid" (faster), "edge", "hausdorff-edge",
#' "least-cost" or "commute-time". If the type is equal to "least-cost" or "commute-time", then you have
#' to use the "resistance" argument. To see more arguments see the ?distancefile.
#' @param metric character. Choose a Betweenness Centrality Metric: "BC" or "BCIIC" considering topologycal distances or "BCPC" considering maximum product probabilities.
#' @param distance_thresholds numeric. Distance or distances thresholds (meters) to establish connections. For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "BPC" metric.
#' @param LA numeric. Maximum landscape attribute (attribute unit, if attribute is NULL then unit is equal to ha).
#' @param coneforpath character. Path to Conefor 2.6 with command line interface
#' (\url{http://www.conefor.org/coneforsensinode.html}). Example, "C:/Users/coneforWin64.exe".
#' @param dA logical. If TRUE, then the delta attribute will be added to the node's importance result.
#' @param dvars logical. If TRUE, then the absolute variation will be added to the node's importance result.
#' @param parallel numeric. Specify the number of cores to use for parallel processing, default = NULL.
#' Parallelize the function using furrr package and multiprocess plan.
#' @param rasterparallel logical. If parallel is FALSE and nodes is a raster then you can use this argument to assign the metrics values to the nodes raster. It is useful when raster resolution is less than 100 m2.
#' @param write character. Write output shapefile, example, "C:/ejemplo.shp".
#' @references Saura, S. and Torne, J. (2012). Conefor 2.6. Universidad Politecnica de Madrid. Available at \url{www.conefor.org}.\cr
#'  Freeman L.C. (1977). Set of Measures of Centrality Based on Betweenness. Sociometry 40: 35-41.\cr
#'  Bodin, O. and Saura, S. (2010). Ranking individual habitat patches as connectivity providers: integrating network analysis and patch removal experiments. Ecological Modelling 221: 2393-2405.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' data("vegetation_patches", package = "Makurhini")
#' nrow(vegetation_patches) #Number of patches
#'
#' #Two distance thresholds
#' BCIIC <- MK_BCentrality(nodes = vegetation_patches, id = "id",
#'             coneforpath = "C:/Users/coneforWin64.exe",
#'             distance = list(type = "centroid"),
#'             metric = "BCIIC", LA = NULL,
#'             distance_thresholds = c(10000, 30000)) #10 and 30 km
#'
#' #Using raster
#' data("raster_vegetation_patches", package = "Makurhini")
#' BCPC <- MK_BCentrality(nodes = raster_vegetation_patches,
#'            coneforpath = "C:/Users/coneforWin64.exe",
#'            attribute = NULL,
#'            distance = list(type = "centroid"),
#'            metric = "BCPC", probability = 0.5,
#'            LA = NULL,
#'            distance_thresholds = 40000) #40 km
#'
#' BCPC_parallel <- MK_BCentrality(nodes = raster_vegetation_patches,
#'                     coneforpath = "C:/Users/coneforWin64.exe",
#'                     id = "id", attribute = NULL,
#'                     distance = list(type = "centroid"),
#'                     metric = "BCPC", LA = NULL, probability = 0.5,
#'                     distance_thresholds = c(40000, 60000),
#'                     parallel = 4) #40 and 60 km
#' }
#' @importFrom sf st_as_sf
#' @importFrom progressr handlers handler_pbcol progressor
#' @importFrom crayon bgWhite white bgCyan
#' @importFrom utils write.table
#' @importFrom raster values as.matrix extent raster stack extent<- writeRaster reclassify crs crs<-
#' @importFrom future multiprocess plan availableCores
#' @importFrom furrr future_map
#' @importFrom purrr walk
#' @importFrom magrittr %>%

MK_BCentrality <- function(nodes, id, attribute  = NULL, area_unit = "ha",
                        distance = list(type= "centroid", resistance = NULL),
                        metric = c("BC", "BCIIC", "BCPC"), distance_thresholds = NULL,
                        probability = NULL, LA = NULL, coneforpath = NULL,
                        dA = FALSE, dvars = FALSE,
                        parallel = FALSE, rasterparallel = FALSE, write = NULL) {
  . = NULL
  if (missing(nodes)) {
    stop("error missing shapefile file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing shapefile file of nodes")
    }
  }

  if (!metric %in% c("BC", "BCIIC", "BCPC")) {
    stop("Type must be either 'BC', 'BCIIC', or 'BCPC'")
  }

  if (isTRUE(unique(metric == c("BC", "BCIIC", "BCPC")))) {
    metric = "BC"
  }

  if (metric == "BCPC") {
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

  if(isFALSE(parallel)){
    parallel <- NULL
  }

  if(isTRUE(parallel)){
    message(paste0("The number of available cores is ", as.numeric(availableCores()),
                   ", so ", as.numeric(availableCores()), " cores will be used."))
    parallel <- as.numeric(availableCores())-2
  }

  if (!is.null(coneforpath)) {
    if (!dir.exists(dirname(coneforpath))) {
      stop("error, output folder does not exist")
    }

    if(!file.exists(coneforpath)){
      stop("error, Conefor 2.6 with command line interface does not exist")
    }
  }

  if (class(nodes)[1] == "SpatialPolygonsDataFrame") {
    nodes <- st_as_sf(nodes)
  }

  options(warn = -1)
  ttt.2 <- getwd(); temp.1 <- paste0(tempdir(), "/TempInputs", sample(1:1000, 1, replace = TRUE))
  dir.create(temp.1, recursive = T)

  if(class(nodes)[1] == "sf"){
    if (is.null(id)) {
      nodes$IdTemp <- 1:nrow(nodes)
      } else {
        nodes$IdTemp <- nodes[[id]]
      }
    id = "IdTemp"
  } else {
    id = NULL
    }


  nodesfile(nodes, id = id, attribute = attribute, area_unit = area_unit,
            write = paste0(temp.1, "/nodes.txt"))


  distancefile(nodes = nodes,  id = id, type = distance$type,
               distance_unit = distance$distance_unit, keep = distance$keep,
               resistance = distance$resistance,
               resist.units = distance$resist.units,
               CostFun = distance$CostFun, ngh = distance$ngh,
               mask = distance$mask,
               threshold = distance$threshold,
               geometry_out = distance$geometry_out,
               bounding_circles = distance$bounding_circles,
               parallel = distance$parallel,
               edgeParallel = distance$edgeParallel, pairwise = TRUE,
               least_cost.java = distance$least_cost.java,
               cores.java = distance$cores.java, ram.java = distance$ram.java,
               write = paste0(temp.1,"/Dist.txt"))

  setwd(temp.1)
  if (is.null(distance$threshold)) {
    pairs = "all"
  } else {
    pairs = "notall"
  }

  x = NULL; loop <- 1:length(distance_thresholds)

  if(length(distance_thresholds)>1){
    handlers(global = T)
    handlers(handler_pbcol(complete = function(s) crayon::bgYellow(crayon::white(s)),
                           incomplete = function(s) crayon::bgWhite(crayon::black(s)),
                           intrusiveness = 2))
    pb <- progressor(along = loop)
  }

  if(isFALSE(parallel)){
    BC_metric <- tryCatch(lapply(loop, function(x){
      x <- distance_thresholds[x]

      if(length(distance_thresholds)>1){
        pb()
      }
        dMetric <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                              coneforpath = coneforpath,
                              typeconnection = "dist",
                              typepairs = pairs, index = metric,
                              thdist = x, multdist = NULL, conprob = probability,
                              onlyoverall = FALSE, LA = LA, nrestauration = FALSE,
                              prefix = NULL, write = NULL)

        if(class(nodes)[1] == "sf"){

          result_interm <- merge_conefor(datat = dMetric[[which(lapply(dMetric, function(x) ncol(x)) >= 11)]],
                                         pattern = NULL,
                                         merge_shape = nodes, id = "IdTemp",
                                         write = if (!is.null(write)) paste0(write, "_d", x,".shp"),
                                         dA = dA, var = dvars)
          result_interm$"IdTemp" <- NULL

        } else {
          datat <- dMetric[[which(lapply(dMetric, function(x) ncol(x)) >= 11)]]
          datat <- as.data.frame(datat);  datat <- datat[, as.numeric(which(colSums(!is.na(datat)) > 0))]
          if(isFALSE(dA)){
            datat$dA <- NULL; datat$varA <- NULL
            datat[,which(unique(is.na(datat)) == TRUE)] <- NULL
          }

          if (isFALSE(dvars)){
            datat <- datat[which(grepl("var",names(datat), fixed = TRUE) == FALSE)]
          }

          ###rasters values
          rp <- unique(raster::values(nodes)) %>% as.vector(.); rp <- rp[which(!is.na(rp))]
          m <- matrix(nrow = nrow(datat), ncol = 2); m[,1] <- datat[,1]
          if(isTRUE(rasterparallel)){
            works <- as.numeric(availableCores())-1
            plan(strategy = multiprocess, gc=TRUE, workers = works)
            r_metric <- tryCatch(future_map(2:ncol(datat), function(c){
              x1 <- datat[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)}, .progress = TRUE), error = function(err) err)
            close_multiprocess(works)
          } else {
            r_metric <- lapply(2:ncol(datat), function(c){
              x1 <- datat[,c(1, c)]
              for(i in rp){
                m[which(m == i),2] <- x1[which(x1[,1]== i),2]
              }
              x1 <- reclassify(nodes, rcl = m)
              return(x1)})
          }

          result_interm <- list()

          for(i in 2:(length(r_metric)+1)){
            result_interm[[i]] <- r_metric[[i-1]]
          }

          result_interm[[1]] <- nodes; result_interm <- stack(result_interm)
          names(result_interm) <- names(datat[,1:ncol(datat)])

          if (!is.null(write)){
            n <- names(datat[,1:ncol(datat)]); walk(as.list(2:length(n)), function(w){
              x1 <- result_interm[[w]]
              crs(x1) <- crs(result_interm)
              writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x, ".tif"), overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
            })
          }
        }
        return(result_interm)
      }), error = function(err) err)
  } else {
    works <- as.numeric(availableCores())-1; works <- if(parallel > works){works}else{parallel}
    plan(strategy = multiprocess, gc = TRUE, workers = works)
    BC_metric <- tryCatch(future_map(loop, function(x) {
      x <- distance_thresholds[x]

      if(length(distance_thresholds)>1){
        pb()
      }
      dMetric <- EstConefor(nodeFile = "nodes.txt", connectionFile = "Dist.txt",
                            coneforpath = coneforpath,
                            typeconnection = "dist", typepairs = pairs, index = metric,
                            thdist = x, multdist = NULL, conprob = probability,
                            onlyoverall = FALSE, LA = LA, nrestauration = FALSE,
                            prefix = NULL, write = NULL)

      if(class(nodes)[1] == "sf"){
        result_interm <- merge_conefor(datat = dMetric[[which(lapply(dMetric, function(x) ncol(x)) >= 11)]], pattern = NULL,
                                       merge_shape = nodes, id = "IdTemp",
                                       write = if (!is.null(write)) paste0(write, "_d", x,".shp"),
                                       dA = dA, var = dvars)
        result_interm$"IdTemp" <- NULL

      } else {
        datat <- dMetric[[which(lapply(dMetric, function(x) ncol(x)) >= 11)]]
        datat <- as.data.frame(datat)
        datat <- datat[, as.numeric(which(colSums(!is.na(datat)) > 0))]
        if(isFALSE(dA)){
          datat$dA <- NULL; datat$varA <- NULL
          datat[,which(unique(is.na(datat)) == TRUE)] <- NULL
        }

        if (isFALSE(dvars)){
          datat <- datat[which(grepl("var",names(datat), fixed = TRUE) == FALSE)]
        }

        ###rasters values
        rp <- unique(raster::values(nodes)) %>% as.vector(.); rp <- rp[which(!is.na(rp))]

        m <- matrix(nrow = nrow(datat), ncol = 2); m[,1] <- datat[,1]
        r_metric <- lapply(2:ncol(datat), function(c){
          x1 <- datat[,c(1, c)]
          for(i in rp){
            m[which(m == i),2] <- x1[which(x1[,1]== i),2]
          }
          x1 <- reclassify(nodes, rcl = m)
          return(x1)})

        result_interm <- list()

        for(i in 2:(length(r_metric)+1)){
          result_interm[[i]] <- r_metric[[i-1]]
        }

        result_interm[[1]] <- nodes; result_interm <- stack(result_interm)
        names(result_interm) <- names(datat[,1:ncol(datat)])

        if (!is.null(write)){
          n <- names(datat[,1:ncol(datat)]); walk(as.list(2:length(n)), function(w){
            x1 <- result_interm[[w]]
            crs(x1) <- crs(result_interm)
            writeRaster(x1, filename = paste0(write, "_", n[w], "_",  x, ".tif"), overwrite = TRUE, options = c("COMPRESS=LZW", "TFW=YES"))
          })
        }
      }
      return(result_interm)
    }), error = function(err) err)
    close_multiprocess(works)
    }

  if (inherits(BC_metric, "error"))  {
    setwd(ttt.2)
    stop(BC_metric)
    } else {
      if (length(distance_thresholds) == 1) {
        BC_metric <- BC_metric[[1]]
        } else {
          names(BC_metric) <- paste0("d", distance_thresholds)
          }
      setwd(ttt.2)
      }
  return(BC_metric)
  }
