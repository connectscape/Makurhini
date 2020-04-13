#' Protected Connected (ProtConn)
#'
#' Use the CONEFOR command line to estimate Protected Connected (ProtConn) indexes.
#' @param nodes object of class sf, sfc, sfg or SpatialPolygons. Example, Protected areas shapefile.
#' @param region object of class sf, sfc, sfg or SpatialPolygons. Region shapefile.
#' @param thintersect numeric. Threshold of intersection in percentage allowed to select or not a target geometry.
#'  Example, if thintersect is equal to 90 then a node will be selected only if the intersection between the node and
#'   the region is >= 90 percentage. If NULL, thintersect will be 0 (default)
#' @param attribute character. Select the nodes attribute: "Intersected area" = Intersected Protected areas; or
#' another specific column name with the nodes attribute, ideally this attribute mus be an area-weighted index,
#'  otherwise the interpretation of the protconn index may change.
#' @param area_unit character. Attribute area units. You can set an area unit, "Makurhini::unit_covert()" compatible unit ("m2", "Dam2, "km2", "ha", "inch2", "foot2", "yard2", "mile2"). Default equal to hectares "ha".
#' @param res_attribute numeric. If the attribute is no equal to "Area" or "Intersected area" then nodes will be converted to raster to extract values in one  process step, you can set the raster resolution, default = 150.
#' @param distance list. See distancefile(). E.g.: list(type= "centroid", resistance = NULL).
#' @param probability numeric. Connection probability to the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection. Use in case of selecting the "PC" metric.
#' @param distance_thresholds numeric. Distance or distances thresholds to establish connections (meters). For example, one distance: distance_threshold = 30000; two or more specific distances:
#'  distance_thresholds = c(30000, 50000); sequence distances: distance_thresholds = seq(10000,100000, 10000).
#' @param transboundary numeric. Buffer to select polygones in a second round, their attribute value = 0, see Saura et al. 2017. You can set one transboundary value or one per each threshold distance.
#' @param LA numeric. Maximum Landscape Attribute.
#' @param keep numeric. Simplification of the region's geometry to speed the select by location ("Makurhini::MK_selectbyloc").
#'  The value can range from 0 to 1 and is the proportion of points to retain (default 1). The lower the value,
#'   the higher the speed but the greater uncertainty.
#' @param plot logical. Plot the main ProtConn indicators and fractions, default = FALSE.
#' @param dPC logical. if TRUE dPC indexes are merged with the nodes shapefile, default = FALSE.
#' @param write character. Output folder including the output file name without extension, e.g., "C:/ProtConn/Protfiles".
#' @param SAGA Logical. Optimize the large process using SAGA GIS and RSAGA package (see, \url{https://github.com/r-spatial/RSAGA}).
#' @param intern logical. Show the progress of the process, default = TRUE.
#' @return
#' Table with the following ProtConn values: ECA, Prot, ProtConn, ProtUnconn, RelConn, ProtConn[design], ProtConn[bound], ProtConn[Prot], ProtConn[Within],
#'  ProtConn[Contig], ProtConn[Trans], ProtConn[Unprot], ProtConn[Within][land], ProtConn[Contig][land],
#'  ProtConn[Unprot][land], ProtConn[Trans][land] \cr
#' \cr
#' *If plot is not NULL a list is returned with the ProtConn table and a plots.
#' @references
#' Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the world’s ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
#' Saura, S., Bertzky, B., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2018). Protected area connectivity: Shortfalls in global targets and country-level priorities. Biological Conservation, 219(October 2017), 53–67.
#' @export
#' @examples
#' \dontrun{
#' library(Makurhini)
#' library(raster)
#' data("Protected_areas", package = "Makurhini")
#' #plot(Protected_areas, col="green")
#'
#' data("regions", package = "Makurhini")
#' region <- regions[1,]
#' plot(region, col="blue")
#'
#' test <- MK_ProtConn(nodes = Protected_areas, region = region,
#'                     attribute = "Intersected area", area_unit = "ha",
#'                     distance = list(type= "centroid"),
#'                     distance_thresholds = c(50000, 10000),
#'                     probability = 0.5, transboundary = 50000,
#'                     LA = NULL, plot = TRUE, dPC = FALSE,
#'                     write = NULL, SAGA = FALSE, intern = TRUE)
#' test$d50000
#' test$d10000
#' }
#' @importFrom sf st_as_sf st_cast st_zm st_simplify st_buffer write_sf st_intersection st_difference st_area
#' @importFrom magrittr %>%
#' @importFrom rmapshaper ms_dissolve ms_simplify
#' @importFrom dplyr progress_estimated mutate group_by summarize
#' @importFrom raster raster extract
#' @importFrom fasterize fasterize
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot geom_bar aes position_dodge labs rel theme_bw theme element_blank element_text scale_fill_manual geom_hline scale_linetype_manual guide_legend margin
#' @importFrom purrr compact
#' @importFrom ggpubr ggarrange
#' @importFrom formattable formattable formatter style color_tile as.htmlwidget
#' @importFrom methods as
#' @importFrom rlang .data
MK_ProtConn <- function(nodes, region, thintersect = NULL,
                        attribute = "Intersected area",
                        area_unit = "ha",
                        res_attribute = 150,
                        distance = list(type= "centroid", resistance = NULL),
                        distance_thresholds, probability,
                        transboundary = NULL, LA = NULL,
                        keep = 1,
                        plot = FALSE,
                        dPC = FALSE, write = NULL, SAGA = FALSE,
                        intern = TRUE){
  options(warn = -1)
  err <- tryCatch(detach("package:plyr", unload = TRUE), error = function(err)err)

  if (missing(nodes)) {
    stop("error missing file of nodes")
  } else {
    if (is.numeric(nodes) | is.character(nodes)) {
      stop("error missing file of nodes")
    }
  }

  if (missing(region)) {
    stop("error missing file of region")
  } else {
    if (is.numeric(region) | is.character(region)) {
      stop("error missing file of region")
    }
  }

  if (is.null(probability) | !is.numeric(probability)) {
    stop("error missing probability")
  }

  if (is.null(distance_thresholds)) {
    stop("error missing numeric distance threshold(s)")
  }

  if (!is.null(write)) {
    if (!dir.exists(dirname(write))) {
      stop("error, output folder does not exist")
    }
  }

  if(is.null(attribute)){
    attribute = "Intersected area"
  }

  if(is.null(thintersect)){
    thintersect = 0
  }

  if(is.character(thintersect)){
    stop("thintersect must be NULL or numeric")
  }

  #Nodes and region topol. correction
  region <- st_as_sf(region) %>% st_zm() %>%
    st_buffer(., dist = 0) %>% st_cast("POLYGON")

  nodes <- st_as_sf(nodes) %>% st_zm() %>%
    st_buffer(., dist = 0) %>% st_cast("POLYGON")

  #1 First selection to reduce the extent of work and dissolve polygones
  select_distance <- max(c(transboundary, distance_thresholds)) * 2

  if(num_vert(region) > 100){
    region_1 <- tryCatch(st_simplify(x = region, dTolerance = 10000, preserveTopology = TRUE),
                         error = function(err)err)

    if (inherits(region_1, "error")){
      region_1 <- st_buffer(x = region, dist = select_distance) %>%
        st_buffer(., dist = 0) %>% st_cast("POLYGON")
    } else {
      region_1 <- st_buffer(x = region_1, dist = select_distance)
    }

  } else {
    region_1 <- st_buffer(x = region, dist = select_distance)
  }

  select_PA <- tryCatch(over_poly(x = nodes, y = region_1), error = function(err)err)

  if (inherits(select_PA, "error")){
    nodes <- st_buffer(x = nodes, dist = 0)
    select_PA <- over_poly(x = nodes, y = region_1)
    nodes.1 <- nodes[!is.na(select_PA),]
  } else {
    nodes.1 <- nodes[!is.na(select_PA),]
  }

  ## Loop transboundary
  if(length(transboundary) > 1){
    pb <- progress_estimated(length(transboundary), 0)
  }

  result_protconn <- tryCatch(map(as.list(if(is.null(transboundary)){100}else{transboundary}), function(x){
    if(nrow(nodes.1) > 1){
      nodes.1$IDTemp <- 1:nrow(nodes.1)
      nodes.1$dis <- 1

      if(num_vert(region) > 100){
        region_2 <- ms_simplify(input = region, keep = keep, keep_shapes = TRUE, explode = TRUE) %>%
          ms_dissolve()
      } else {
        region_2 <- region
      }

      nodes.2 <- MK_selectbyloc(target = nodes.1, sourcelyr = region_2,
                                selreg = "M2", thintersect = thintersect,
                                area_unit = area_unit,
                                transboundary = if(is.null(x)){100}else{x},
                                SAGA = SAGA)

      if(nrow(nodes.2) > 0){
        nodes.2 <- ms_dissolve(nodes.2, field = "dis") %>% st_buffer(., dist = 0) %>%
          st_cast("POLYGON")
        nodes.2 <- MK_selectbyloc(target = nodes.2, sourcelyr = region_2,
                                  selreg = "M2", thintersect = thintersect,
                                  area_unit = area_unit,
                                  transboundary = if(is.null(x)){100}else{x},
                                  SAGA = SAGA)
      }

      nodes.2$dis <- NULL

      if(is.null(x)){
        nodes.2 <- nodes.2[nodes.2$transboundary == 1,]
      }

      if(nrow(nodes.2) > 0){
        if(!is.null(x) & attribute == "Intersected area"){
          nodes.2t1 <- st_intersection(nodes.2[nodes.2$transboundary==1,], region_2) %>%
            st_cast("POLYGON")
          nodes.2t1$rmapshaperid <- NULL
          nodes.2t1$TempID <- NULL
          nodes.2t1$Area2 <- unit_convert(as.numeric(st_area(nodes.2t1)), "m2", area_unit)

          nodes.2t2 <- st_difference(nodes.2[nodes.2$transboundary==1,], region_2) %>%
            st_cast("POLYGON")
          if(nrow(nodes.2t2) > 0){
            nodes.2t2$rmapshaperid <- NULL
            nodes.2t2$TempID <- NULL
            nodes.2t2$transboundary <- 0
            nodes.2t2 <- rbind(nodes.2[nodes.2$transboundary==0,], nodes.2t2)
            nodes.2 <- rbind(nodes.2t1, nodes.2t2)
          }
        }
      }
      area2_sort <- sort(x = nodes.2[nodes.2$transboundary == 1,]$Area2, decreasing = T)
    } else if (nrow(nodes.1) == 1){
      if(nrow(nodes.1) == 1){
        nodes.2 <- nodes.1
        nodes.2$transboundary <- 1
        nodes.2$Area2 <- unit_convert(as.numeric(st_area(nodes.2)), "m2", area_unit)
        area2_sort <- nodes.2$Area2
      } else {
        nodes.2 <- nodes.1
      }
    } else {
      nodes.2 <- nodes.1
    }

    #Step 2. Distances
    #If there are protected areas get distance file
    if(nrow(nodes.2) > 1 ){
      if(area2_sort[[1]] < sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit))){
        nodes.2$IDTemp <- 1:nrow(nodes.2)
        bound_nodes <- nodes.2

        ##Nodes region
        n <- seq(1, nrow(region))
        bound_nodes2 <- region
        bound_nodes2$IDTemp <- n + nrow(bound_nodes) #IDs
        bound_nodes2 <- bound_nodes2[,"IDTemp"]

        #colocar los mismos campos
        bound_nodes2$rmapshaperid <- 0
        bound_nodes2$Area1 <- 0
        bound_nodes2$Area2 <- 0
        bound_nodes2$PercPi <- 0
        bound_nodes2$transboundary <- 0
        bound_nodes2$AreaTemp <- 0

        bound_nodes <- rbind(bound_nodes, bound_nodes2[, names(bound_nodes)])
        bound_nodes[which(is.na(bound_nodes$rmapshaperid)), names(bound_nodes)] <- 0

        #If cost or commutate distance will be used
        if(!is.null(distance$resistance) & is.null(distance$mask)){
          mask <- region_1
          resistance_protconn <- crop(distance$resistance, as(st_buffer(mask, dist = distance_thresholds), 'Spatial'))
          distance$resistance <- resistance_protconn
        }

        distance_base <- tryCatch(distancefile(nodes = bound_nodes,  id = "IDTemp", type = distance$type,
                                               distance_unit = distance$distance_unit, keep = distance$keep,
                                               resistance = distance$resistance,
                                               CostFun = distance$CostFun, ngh = distance$ngh,
                                               mask = distance$mask,
                                               threshold = distance$threshold,
                                               geometry_out = distance$geometry_out,
                                               bounding_circles = distance$bounding_circles,
                                               parallel = distance$parallel,
                                               edgeParallel = distance$edgeParallel, pairwise = FALSE,
                                               write = NULL), error = function(err)err)

        if (inherits(distance_base, "error")){
          bound_nodes <- st_buffer(bound_nodes, dist = 0)
          distance_base <- tryCatch(distancefile(nodes = bound_nodes,  id = "IDTemp", type = distance$type,
                                                 distance_unit = distance$distance_unit, keep = distance$keep,
                                                 resistance = distance$resistance,
                                                 CostFun = distance$CostFun, ngh = distance$ngh,
                                                 mask = distance$mask,
                                                 threshold = distance$threshold,
                                                 geometry_out = distance$geometry_out,
                                                 bounding_circles = distance$bounding_circles,
                                                 parallel = distance$parallel,
                                                 edgeParallel = distance$edgeParallel, pairwise = FALSE,
                                                 write = NULL), error = function(err)err)
          if (inherits(distance_base, "error")){
            stop("distance file error")
          }
        }

        #Correction of centroid bound
        if(length(bound_nodes2$IDTemp) == 1){
          distance_base[which(row.names(distance_base) == as.numeric(bound_nodes2$IDTemp)),] <- 0
          distance_base[,which(colnames(distance_base) == as.numeric(bound_nodes2$IDTemp))] <- 0
        } else {
          for(i in bound_nodes2$IDTemp){
            over_nodes <- tryCatch(over_poly(x = nodes.2, y = bound_nodes2[which(bound_nodes2$IDTemp == i),]), error = function(err)err)
            over_nodes <- nodes.2[!is.na(over_nodes),]
            distance_base[which(row.names(distance_base) == i),
                          which(row.names(distance_base) %in% as.numeric(over_nodes$IDTemp))] <-0
            distance_base[which(colnames(distance_base) %in% as.numeric(over_nodes$IDTemp)),
                          which(row.names(distance_base) == i)] <-0
          }
        }
      }
    }


    nodes_base <- nodes.2

    #Step 3. Loop distance thresholds
    #progress
    i = NULL
    pb2 <- progress_estimated(length(distance_thresholds), 0)

    result_protconn_2 <- tryCatch(map(as.list(distance_thresholds), function(i){
      nodes.2 <- nodes_base
      #If there are protected areas
      if (nrow(nodes.2) > 1 ){
        if(area2_sort[[1]] < sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit))){
          if (attribute == "Intersected area"){
            nodes.2$AreaTemp <- nodes.2$Area2 * nodes.2$transboundary
          } else if (attribute == "Area"){
            nodes.2$AreaTemp <- nodes.2$Area1 * nodes.2$transboundary
          } else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
            r.raster <- raster(nodes.1, res= res_attribute)
            r.raster <- fasterize(nodes.1, r.raster, field = attribute, background = 0)
            Mode <- function(x, na.rm = FALSE) {
              if(na.rm){
                x = subset(x, !is.na(x))
              }
              ux <- unique(x)
              ux <- ux[which.max(tabulate(match(x, ux)))]
              if(is.na(ux)){
                ux = 0
              }
              return(ux)
            }
            majority <- extract(r.raster, nodes.2, fun = Mode)
            majority[is.na(majority)] <- 0
            nodes.2$AreaTemp <- majority * nodes.2$transboundary
          } else {
            stop("missing attribute")
          }

          #Get ECA with Transboundary PAs
          nodes_protconn <- tryCatch(nodesfile(nodes = nodes.2, id = "IDTemp", attribute = "AreaTemp",
                                               write = NULL), error = function(err)err)


          if (inherits(nodes_protconn, "error")){
            stop("nodes file error")
          }

          `%nin%` = Negate(`%in%`)
          distance_protconn <- distance_base[which(row.names(distance_base) %nin% as.numeric(bound_nodes2$IDTemp)),
                                             which(colnames(distance_base) %nin% as.numeric(bound_nodes2$IDTemp))]

          #Landscape attribute for ECA
          if (is.null(LA) & isTRUE(attribute %in% c("Intersected area"))){
            LA <- sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit))
          } else {
            if(is.null(LA) & isFALSE(attribute %in% c("Area", "Intersected area"))){
              stop("misssing LA, LA is necessary when you choose attribute a different to Area and Intersected area")
            } else {
              LA <- LA
            }
          }

          #Conefor command line
          result.1 <- tryCatch(MK_dPCIIC(nodes = nodes_protconn, attribute = names(nodes_protconn)[2],
                                         distance = distance_protconn,
                                         metric = "PC", probability = probability, LA = LA,
                                         overall = TRUE, onlyoverall = if(isTRUE(dPC)){FALSE}else{TRUE},
                                         distance_thresholds = i), error = function(err)err)

          if (inherits(result.1, "error")){
            stop(result.1)
          }

          #dPC merge
          if(isTRUE(dPC)){
            rdpc <- result.1[[1]]
            if(length(unique(rdpc[,4])) > 1 | sum(rdpc[,4]) != 0){
              dPC.2 <- base::merge(nodes.2, rdpc[,c( which(names(rdpc) == "IDTemp"),
                                                     which(names(rdpc) == "dPC"):ncol(rdpc))], by="IDTemp")
              dPC.2$"IDTemp" <- NULL
              dPC.2$OBJECTID <- 1:nrow(dPC.2)
              dPC.2 <- dPC.2[, c("OBJECTID", "PercPi", "transboundary", "AreaTemp", "dPC", "dPCintra", "dPCflux", "dPCconnector")]
            } else if (unique(rdpc[,4]) == 0 ){
              dPC.2 <- "dPC index equal 0. Possible explanation: Only transboundary protected areas where selected, please, check the thintersect parameter."
            } else {
              dPC.2 <- "dPC index equal 0. Possible explanation: Only transboundary protected areas where selected, please, check the thintersect parameter."
            }
          }

          #Initial ProtConn table
          #Landscape attribute for protected land
          data <- if(isTRUE(dPC)){result.1[[2]]} else {result.1}

          DataProtconn <- as.data.frame(cbind(ECA = data[[2,2]], PC = data[[3,2]], LA))
          #
          if (attribute == "Intersected area"){
            DataProtconn$a <- sum(nodes.2$Area2 * nodes.2$transboundary)
          } else if (attribute == "Area"){
            DataProtconn$a <- sum(nodes.2$Area1 * nodes.2$transboundary)
          } else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
            DataProtconn$a <- sum(nodes.2$Area2 * nodes.2$transboundary)
          }

          #ProtConn Indicators
          DataProtconn$Prot <- 100 * (DataProtconn$a / sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit)))#Percentage of the study region covered by PAs
          DataProtconn$ProtConn <- 100 * (DataProtconn$ECA / DataProtconn$LA)
          DataProtconn$ProtUnconn <- DataProtconn$Prot - DataProtconn$ProtConn
          DataProtconn$RelConn <- 100 * (DataProtconn$ProtConn / DataProtconn$Prot)
          DataProtconn$Unprotected <- 100 - DataProtconn$Prot

          #Get ECA WITHOUT Transboundary PAs
          #If not Transboundary PAs. 1 = No Transboundary Protected Area
          if(isFALSE(unique(nodes.2$transboundary) == 1)){
            ECA2 <- as.data.frame(cbind(ECA = 0, PC = 0))
            AreaTemp1 <- 0
            AreaTemp2 <- 0
          } else {
            #If Transboundary PAs
            nodes.2 <- nodes.2[nodes.2$transboundary == 1,]

            if (nrow(nodes.2) > 1){
              #Nodes file
              nodes_protconn.2 <-  nodesfile(nodes = nodes.2, id = "IDTemp", attribute = "AreaTemp",
                                             write = NULL)

              #Distance file
              distance_protconn.2 <- distance_base[which(row.names(distance_base) %in% as.numeric(nodes.2$IDTemp)),
                                                   which(colnames(distance_base) %in% as.numeric(nodes.2$IDTemp))]

              result.2 <- tryCatch(MK_dPCIIC(nodes = nodes_protconn.2, attribute = names(nodes_protconn.2)[2],
                                             distance = distance_protconn.2,
                                             metric = "PC", probability = probability, LA = LA,
                                             overall = TRUE, onlyoverall = TRUE,
                                             distance_thresholds = i), error = function(err)err)


              if (inherits(result.2, "error")){
                stop(result.2)
              }

              ECA2 <- as.data.frame(cbind(ECA = result.2[[2,2]], PC = result.2[[3,2]]))

              #a with aggregated polygones and exploded and withouth transboundary attribute = 0
              AreaTemp1 <- nodes.2$AreaTemp
              #Get a'(not dissolved polygones)
              TopologError <- tryCatch(over_poly(nodes.1, nodes.2), error = function(err)err)
              nodes.3 <- nodes.1[!is.na(TopologError),]

              if (inherits(TopologError, "error")){
                nodes.1 <- st_buffer(nodes.1, dist = 0) %>% st_cast("POLYGON")
                nodes.2 <- st_buffer(nodes.2, dist = 0) %>% st_cast("POLYGON")
                TopologError <- tryCatch(over_poly(nodes.1, nodes.2), error = function(err)err)
                nodes.3 <- nodes.1[!is.na(TopologError),]
              }
              nodes.3 <- nodes.3[,"IDTemp"]
              nodes.3$IDTemp <- 1:nrow(nodes.3)

              if(nrow(nodes.3) >= 1){
                nodes.3$A <- unit_convert(as.numeric(st_area(nodes.3, by_element = TRUE)), "m2", area_unit)
                pi2 <- tryCatch(st_intersection(nodes.3, region_2), error = function(err)err)

                if (inherits(pi2, "error")){
                  nodes.3 <- st_buffer(nodes.3, dist = 0)
                  region_2 <-  st_buffer(region_2, dist = 0)
                  pi2 <- st_intersection(nodes.3, region_2)
                }

                attArea2 <- pi2 %>%
                  mutate(A2 = unit_convert(as.numeric(st_area(.)), "m2", area_unit) %>% as.numeric())%>%
                  as_tibble() %>% group_by(.data$IDTemp) %>%
                  dplyr::summarize(Area1 = sum(.data$A), Area2 = sum(.data$A2))

                attArea2$PercPi <- as.numeric((attArea2$Area2 * 100) / attArea2$Area1)
                nodes.3 <- base::merge(nodes.3, attArea2, by = "IDTemp", all = TRUE)

                if(length(which(is.na(nodes.3$PercPi))) >= 1){
                  nodes.3$PercPi[is.na(nodes.3$PercPi)] <- 0
                  nodes.3$Area2[is.na(nodes.3$Area2)] <- 0
                  nodes.3$Area1 <- nodes.3$A
                }

                if(thintersect > 0){
                  nodes.3 <- nodes.3[nodes.3$PercPi >= thintersect,]
                }

                if (attribute == "Intersected area"){
                  AreaTemp2 <- nodes.3$Area2
                } else if (attribute =="Area"){
                  AreaTemp2 <- nodes.3$A
                }else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
                  majority <- r.raster$extract(sp=nodes.3, fun = Mode)
                  majority[is.na(majority)] <- 0
                  AreaTemp2 <- majority
                } else {
                  stop("missing attribute")
                }

              } else {
                AreaTemp2 <- 0
              }

            } else if (nrow(nodes.2) == 1){
              ECA2 <- as.data.frame(cbind(ECA = 0, PC = 0))
              AreaTemp1 <- nodes.2$AreaTemp
              #Get a'los poligonos sin agregar
              TopologError <- tryCatch(over_poly(nodes.1, nodes.2), error = function(err)err)
              nodes.3 <- nodes.1[!is.na(TopologError),]

              if (inherits(TopologError, "error")){
                nodes.1 <- st_buffer(nodes.1, dist = 0) %>% st_cast("POLYGON")
                nodes.2 <- st_buffer(nodes.2, dist = 0) %>% st_cast("POLYGON")
                TopologError <- tryCatch(over_poly(nodes.1, nodes.2), error = function(err)err)
                nodes.3 <- nodes.1[!is.na(TopologError),]
              }
              nodes.3 <- nodes.3[,"IDTemp"]
              nodes.3$IDTemp <- 1:nrow(nodes.3)
              if(nrow(nodes.3) >= 1){
                nodes.3$A <- unit_convert(as.numeric(st_area(nodes.3, by_element = TRUE)), "m2", area_unit)
                pi2 <- tryCatch(st_intersection(nodes.3, region_2), error = function(err)err)

                if (inherits(pi2, "error")){
                  nodes.3 <- st_buffer(nodes.3, dist = 0)
                  region_2 <-  st_buffer(region_2, dist = 0)
                  pi2 <- st_intersection(nodes.3, region_2)
                }

                attArea2 <- pi2 %>%
                  mutate(A2 = unit_convert(as.numeric(st_area(.)), "m2", area_unit) %>% as.numeric())%>%
                  as_tibble() %>% group_by(.data$IDTemp) %>%
                  dplyr::summarize(Area1 = sum(.data$A), Area2 = sum(.data$A2))

                attArea2$PercPi <- as.numeric((attArea2$Area2 * 100) / attArea2$Area1)
                nodes.3 <- base::merge(nodes.3, attArea2, by = "IDTemp", all = TRUE)

                if(length(which(is.na(nodes.3$PercPi))) >= 1){
                  nodes.3$PercPi[is.na(nodes.3$PercPi)] <- 0
                  nodes.3$Area2[is.na(nodes.3$Area2)] <- 0
                  nodes.3$Area1 <- nodes.3$A
                }

                if(thintersect > 0){
                  nodes.3 <- nodes.3[nodes.3$PercPi >= thintersect,]
                }

                if (attribute == "Intersected area"){
                  AreaTemp2 <- nodes.3$Area2
                } else if (attribute =="Area"){
                  AreaTemp2 <- nodes.3$A
                } else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
                  majority <- r.raster$extract(sp=nodes.3, fun = Mode)
                  majority[is.na(majority)] <- 0
                  AreaTemp2 <- majority
                } else {
                  stop("missing attribute")
                }

              } else {
                AreaTemp2 = 0
              }

            } else{
              ECA2 <- as.data.frame(cbind(ECA = 0, PC = 0))
              AreaTemp1 = 0
              AreaTemp2 = 0
            }
          }

          ###Proconn bound
          if (nrow(bound_nodes) > 1){
            if (attribute == "Intersected area"){
              bound_nodes$AreaTemp <- bound_nodes$Area2 * bound_nodes$transboundary
            } else if (attribute == "Area"){
              bound_nodes$AreaTemp <- bound_nodes$Area1 * bound_nodes$transboundary
            }else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
              majority <- r.raster$extract(sp=bound_nodes, fun = Mode)
              majority[is.na(majority)] <- 0
              bound_nodes$AreaTemp <- majority * bound_nodes$transboundary
            } else {
              stop("missing attribute")
            }

            #Nodes file
            nodes_protconn.3 <-  nodesfile(nodes = bound_nodes, id = "IDTemp", attribute = "AreaTemp", write = NULL)

            #Distance file
            result.3 <- tryCatch(MK_dPCIIC(nodes = nodes_protconn.3, attribute = names(nodes_protconn.3)[2],
                                           distance = distance_base,
                                           metric = "PC", probability = probability, LA = LA,
                                           overall = TRUE, onlyoverall = TRUE,
                                           distance_thresholds = i), error = function(err)err)

            if (inherits(result.3, "error")){
              stop(result.3)
            }
            ECA3 <- as.data.frame(cbind(ECA = result.3[[2,2]], PC = result.3[[3,2]]))

            DataProtconn$ProConn_design <- (100*(ECA3$ECA/DataProtconn$LA)) - (100 * (DataProtconn$ECA / DataProtconn$LA))
            DataProtconn$ProConn_Bound<- DataProtconn$Prot - DataProtconn$ProConn_design
          }

          ###ProtConn fractions
          DataProtconn[is.na(DataProtconn)] <- 0
          ProtConn <- DataProtconn$ProtConn
          DataProtconn$ProtConn_Prot <- ((sqrt(sum(AreaTemp1^2)) / LA) * 100) / ProtConn * 100
          DataProtconn$ProtConn_Trans <- 100 * ((100 * ((DataProtconn$ECA - ECA2$ECA) / LA)) / ProtConn)
          DataProtconn[is.na(DataProtconn)] <- 0

          if ((100 * ((DataProtconn$ECA - ECA2$ECA) / LA)) == ProtConn){
            DataProtconn$ProtConn_Trans <- 0
          }

          if (DataProtconn$ProtConn_Trans == 100){
            DataProtconn$ProtConn_Trans <- 0
          }

          DataProtconn$ProtConn_Unprot <- 100 - DataProtconn$ProtConn_Prot - DataProtconn$ProtConn_Trans

          ########ProtConn[Prot] fractions
          within1 <- ((sqrt(sum(AreaTemp2^2))) / LA) * 100
          within2 <- sqrt(sum(AreaTemp1) / sum(AreaTemp2))
          within3 <- within1 * within2

          DataProtconn$ProtConn_Within <- 100 * (within3/ProtConn)

          DataProtconn$ProtConn_Contig <- DataProtconn$ProtConn_Prot - DataProtconn$ProtConn_Within


          ########As a percentage of the total land area in the study region
          DataProtconn$ProtConn_Within_land <- DataProtconn$ProtConn_Within * (ProtConn / 100)
          DataProtconn$ProtConn_Contig_land <- DataProtconn$ProtConn_Contig * (ProtConn / 100)
          DataProtconn$ProtConn_Unprot_land <- DataProtconn$ProtConn_Unprot * (ProtConn / 100)
          DataProtconn$ProtConn_Trans_land <- DataProtconn$ProtConn_Trans * (ProtConn / 100)

        } else {
          #PA equal region area
          if (is.null(LA) & isTRUE(attribute %in% c("Area", "Intersected area"))){
            LA <- sum(unit_convert(as.numeric(st_area(region))), "m2", area_unit)
          } else {
            if(is.null(LA) & isFALSE(attribute %in% c("Area", "Intersected area"))){
              stop("misssing LA, LA is necessary when you choose attribute a different to Area and Intersected area")
            } else {
              LA <- LA
            }
          }

          if (attribute == "Intersected area"){
            Area <- nodes.2$Area2 * nodes.2$transboundary
          } else if (attribute == "Area"){
            Area <- nodes.2$Area1 * nodes.2$transboundary
          } else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
            r.raster <- raster(nodes.1, res= res_attribute)
            r.raster <- fasterize(nodes.1, r.raster, field = attribute, background = 0)
            majority <- extract(r.raster, nodes.2, fun = Mode)
            Mode <- function(x, na.rm = FALSE) {
              if(na.rm){
                x = subset(x, !is.na(x))
              }
              ux <- unique(x)
              ux <- ux[which.max(tabulate(match(x, ux)))]
              if(is.na(ux)){
                ux = 0
              }
              return(ux)
            }

            majority[is.na(majority)] <- 0
            Area <- majority * nodes.2$transboundary
          } else {
            stop("missing attribute")
          }


          DataProtconn <- data.frame(ECA = if(sum(area2_sort)>LA){LA}else{sum(area2_sort)}, PC = NA, LA = LA, a = if(sum(area2_sort)>LA){LA}else{sum(area2_sort)},
                                     Prot = if(round(100 * (sum(Area) / LA), 3) > 100){100}else{round(100 * (sum(Area) / LA), 3)},
                                     ProtConn = if(round(100 * (sum(Area) / LA), 3) > 100){100}else{round(100 * (sum(Area) / LA), 3)},
                                     ProtUnconn = NA,
                                     RelConn = NA,
                                     Unprotected = if(round(100 - (100 * (sum(Area) / LA)),3) < 0){0}else{round(100 - (100 * (sum(Area) / LA)),3)},
                                     ProConn_design = if(round(100 * (sum(Area) / LA), 3) > 100){100}else{round(100 * (sum(Area) / LA), 3)},
                                     ProConn_Bound = if(round(100 * (sum(Area) / LA), 3) > 100){100}else{round(100 * (sum(Area) / LA), 3)},
                                     ProtConn_Prot = 100,
                                     ProtConn_Trans = NA,
                                     ProtConn_Unprot = NA,
                                     ProtConn_Within = 100,
                                     ProtConn_Contig = NA,
                                     ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                     ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)

          if(isTRUE(dPC)){
            dPC.2 <- "dPC index equal 0. Possible explanation: No protected area was selected."
          }
        }

        #If only one Protected area was selected
      } else if (nrow(nodes.2) == 1){

        if (is.null(LA) & isTRUE(attribute %in% c("Area", "Intersected area"))){
          LA <- sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit))
        } else {
          if(is.null(LA) & isFALSE(attribute %in% c("Area", "Intersected area"))){
            stop("misssing LA, LA is necessary when you choose attribute a different to Area and Intersected area")
          } else {
            LA <- LA
          }
        }
        if (attribute == "Intersected area"){
          Area <- nodes.2$Area2 * nodes.2$transboundary
        } else if (attribute == "Area"){
          Area <- nodes.2$Area1 * nodes.2$transboundary
        } else if (is.character(attribute) & isFALSE(attribute %in% c("Area", "Intersected area"))){
          r.raster <- raster(nodes.1, res= res_attribute)
          r.raster <- fasterize(nodes.1, r.raster, field = attribute, background = 0)
          Mode <- function(x, na.rm = FALSE) {
            if(na.rm){
              x = subset(x, !is.na(x))
            }
            ux <- unique(x)
            ux <- ux[which.max(tabulate(match(x, ux)))]
            if(is.na(ux)){
              ux = 0
            }
            return(ux)
          }
          majority <- extract(r.raster, nodes.2, fun = Mode)
          majority[is.na(majority)] <- 0
          Area <- majority * nodes.2$transboundary
        } else {
          stop("missing attribute")
        }

        DataProtconn <- data.frame(ECA = if(sum(Area)> LA){LA}else{sum(Area)}, PC = NA, LA = LA, a = if(sum(Area)> LA){LA}else{sum(Area)},
                                   Prot = if(round(100 * (sum(Area) / LA), 3)>100){100}else{round(100 * (sum(Area) / LA), 3)},
                                   ProtConn = if(round(100 * (sum(Area) / LA), 3)>100){100}else{round(100 * (sum(Area) / LA), 3)},
                                   ProtUnconn = NA,
                                   RelConn = NA,
                                   Unprotected = if(round(100 - (100 * (sum(Area) / LA)),3) < 0){0}else{round(100 - (100 * (sum(Area) / LA)),3)},
                                   ProConn_design = if(round(100 * (sum(Area) / LA), 3)>100){100}else{round(100 * (sum(Area) / LA), 3)},
                                   ProConn_Bound = if(round(100 * (sum(Area) / LA), 3)>100){100}else{round(100 * (sum(Area) / LA), 3)},
                                   ProtConn_Prot = 100,
                                   ProtConn_Trans = NA,
                                   ProtConn_Unprot = NA,
                                   ProtConn_Within = 100,
                                   ProtConn_Contig = NA,
                                   ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                   ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)


        if(isTRUE(dPC)){
          dPC.2 <- "dPC index equal 0. Possible explanation: No protected area was selected."
        }

        ##################################################
        #If no PA was selected
      } else if (nrow(nodes.2) == 0){
        if (is.null(LA) & isTRUE(attribute %in% c("Area", "Intersected area"))){
          LA <- sum(unit_convert(as.numeric(st_area(region)), "m2", area_unit))
        } else {
          if(is.null(LA) & isFALSE(attribute %in% c("Area", "Intersected area"))){
            stop("misssing LA, LA is necessary when you choose attribute a different to Area and Intersected area")
          } else {
            LA <- LA
          }
        }

        DataProtconn <- data.frame(ECA = NA, PC = NA, LA = LA, a = NA,
                                   Prot = NA,
                                   ProtConn = NA,
                                   ProtUnconn = NA,
                                   RelConn = NA,
                                   Unprotected = 100,
                                   ProConn_design = NA,
                                   ProConn_Bound = NA,
                                   ProtConn_Prot = NA,
                                   ProtConn_Trans = NA,
                                   ProtConn_Unprot = NA,
                                   ProtConn_Within = NA,
                                   ProtConn_Contig = NA,
                                   ProtConn_Within_land = NA, ProtConn_Contig_land = NA,
                                   ProtConn_Unprot_land = NA, ProtConn_Trans_land = NA)


        if(isTRUE(dPC)){
          dPC.2 <- "dPC index equal 0. Possible explanation: No protected area was selected."
        }
      }
      names(DataProtconn)[1:2] <-c("ECA", "PC")

      DataProtconn$LA <- DataProtconn$Unprotected
      DataProtconn[c(4,9)] <- NULL
      names(DataProtconn)[3] <- "Unprotected"
      DataProtconn[is.na(DataProtconn)] <- 0
      result <- list()

      #Plot ProtConn
      if (nrow(nodes.2) > 0 & isTRUE(plot)){
        plot_protconn <- list()
        Data_plot <- as.data.frame(t(round(DataProtconn, 3)))
        Data_plot$name <- row.names(Data_plot)
        names(Data_plot) <- c("Value", "Indicator")
        rownames(Data_plot)<- NULL
        Data_plot$Indicator <- as.character(Data_plot$Indicator)
        Data_plot_1 <- Data_plot[which(Data_plot$Indicator %in% c("Unprotected", "Prot","ProtConn")),]
        Data_plot_1$Indicator <-  c("Unprotected", "Protected", "Protected connected")
        Data_plot_1$Indicator <- factor(Data_plot_1$Indicator, levels = c(Data_plot_1$Indicator[1],
                                                                          Data_plot_1$Indicator[2], Data_plot_1$Indicator[3]))
        names(Data_plot_1) <- c("Values", "name")
        Data_plot_1$col <- c("#C34D51", "#53A768", "#4C72AF")
        Data_plot_1[which(Data_plot_1$Values == 0), ] <- NULL

        if(nrow(Data_plot_1) > 1){
          plot_protconn1 <- ggplot(Data_plot_1, aes(x = Data_plot_1$name, y = Data_plot_1$Values,
                                                    fill = Data_plot_1$name)) +
            geom_bar(position = position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
            labs(title = paste0("ProtConn Indicators: ", i), x = "", y = "Percentage (%)", size = rel(1.2)) +
            theme_bw()  +
            theme(plot.title = element_text(color = "#252525", size = rel(1.4), hjust = 0.5, face = "bold"),
                  axis.title = element_text(color = "#252525", size = rel(1.2)),
                  legend.title = element_blank(),
                  legend.text = element_text(colour = "#252525", size = rel(1.2)),
                  axis.text= element_text(colour = "#525252", size = rel(1)))+
            scale_fill_manual(values = Data_plot_1$col) +
            geom_hline(aes(yintercept = 17, linetype = "Aichi Target (17%)"), colour = 'black', size = 1.2) +
            scale_linetype_manual(name = " Aichi Target", values = c(2, 2),
                                  guide = guide_legend(override.aes = list(color = c("black"), size = 0.8)))

          plot_protconn[[1]] <- plot_protconn1
        }

        Data_plot_2 <- Data_plot[which(Data_plot$Indicator %in% c("ProtConn_Trans", "ProtConn_Unprot", "ProtConn_Within", "ProtConn_Contig")),]
        Data_plot_2$Indicator <-  c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]")
        Data_plot_2$Indicator <- factor(Data_plot_2$Indicator, levels = c(Data_plot_2$Indicator[3], Data_plot_2$Indicator[4], Data_plot_2$Indicator[2], Data_plot_2$Indicator[1]))
        names(Data_plot_2) <- c("Values", "name")
        Data_plot_2$col <- c("#253494", "#2c7fb8", "#41b6c4", "#7fcdbb")
        Data_plot_2 <- Data_plot_2[which(Data_plot_2$Values > 0), ]

        if(nrow(Data_plot_2) > 1){
          plot_protconn2 <- ggplot(Data_plot_2, aes(x = Data_plot_2$name, y = Data_plot_2$Values,
                                                    fill = Data_plot_2$name)) +
            geom_bar(position = position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
            labs(title = "Protected connected fraction", x = "", y = "Percentage (%)", size = rel(1.2)) +
            theme_bw()  +
            theme(plot.title = element_text(color = "#252525", size = rel(1.4), hjust = 0.45, face = "bold"),
                  axis.title= element_text(color = "#252525", size = rel(1.2)),
                  plot.margin = margin(0, 5.3, 0, 0.3, "cm"),
                  legend.title= element_blank(),
                  legend.text = element_text(colour = "#252525", size = rel(1.2)),
                  axis.text= element_text(colour = "#525252", size = rel(1)))+
            scale_fill_manual(values = Data_plot_2$col)
          plot_protconn[[2]] <- plot_protconn2
        }
        plot_protconn <- compact(plot_protconn)

        if(length(plot_protconn) == 2){
          figure <- ggarrange(plot_protconn[[1]], plot_protconn[[2]],
                              ncol = 1, nrow = 2)
        } else if (length(plot_protconn) == 1) {
          figure <- plot_protconn[[1]]
        } else {
          figure <- "There are insufficient data to plot"
        }

        result[[2]] <- figure
        names(result)[2] <- "ProtConn Plot"
      }

      if(length(which(DataProtconn < 0)) > 0){
        DataProtconn[,which(DataProtconn < 0)] <- 0
      }
      if(length(which(DataProtconn[3:ncol(DataProtconn)] > 100)) > 0){
        DataProtconn[,which(DataProtconn[3:ncol(DataProtconn)] > 100) + 2] <- 100
      }
      DataProtconn_2 <- t(DataProtconn) %>% as.data.frame()
      DataProtconn_2$Indicator <- row.names(DataProtconn_2)
      DataProtconn_2$Indicator[1] <- "EC(PC)"
      DataProtconn_2$Index <- DataProtconn_2[c(1:2),2]
      DataProtconn_2$Value <- DataProtconn_2[c(1:2),1]

      Value <- DataProtconn_2[c(1:2),1]
      if(Value[1]%%1==0){
        Value <- c(formatC(as.numeric(Value[1]), format="d"),
                   formatC(as.numeric(Value[2]), format="e"))
      } else {
        Value <- c(formatC(as.numeric(Value[1]), format="f", digits = 2),
                   formatC(as.numeric(Value[2]), format="e"))
      }
      DataProtconn_2$Value <- Value
      DataProtconn_2[,1] <- round(DataProtconn_2[,1], 3)
      rownames(DataProtconn_2) <- NULL
      #
      DataProtconn_3 <- DataProtconn_2[3:18,c(3:4, 2, 1)]
      names(DataProtconn_3)[3:4] <- c("ProtConn indicator", "Percentage")
      DataProtconn_3[3:16, 1:2] <- " "
      rownames(DataProtconn_3) <- NULL

      DataProtconn_4 <- formattable(DataProtconn_3, align = c("l","c"),
                                    list(`Index` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `ProtConn indicator` = formatter("span", style = ~ style(color = "#636363", font.weight = "bold")),
                                         `Percentage` = color_tile("white", "orange")))

      result[[1]] <- DataProtconn_4
      names(result)[1] <- "Protected Connected (Viewer Panel)"

      if(isTRUE(dPC) & nrow(nodes.2) > 0){
        result[[3]] <- dPC.2
        names(result)[3] <- "dPC Protected Areas"
      }

      if(!is.null(write)){
        write.csv(result[[1]], paste0(write, "_d", i, "_TableProtConn.csv"), row.names = FALSE)

        if (nrow(nodes.2) > 0 & isTRUE(plot)){
          if(!is.character(figure)){
            tiff(paste0(write, "_d", i, '_ProtConn_plot.tif'), width = 806, height = 641)
            print(figure)
            dev.off()
          } }
        if(isTRUE(dPC)){
          write_sf(result[[3]], paste0(write,"_d", i, "_dPC.shp"), delete_layer = TRUE)
        }
      }
      result2 <- compact(result)

      if (isTRUE(intern) & length(distance_thresholds)>1){
        pb2$tick()$print()
      }
      return(result2)
    }), error = function(err)err)

    #Step 4
    if (isTRUE(intern) & length(transboundary)>1){
      pb$tick()$print()
    }

    if (inherits(result_protconn_2, "error")){
      stop(result_protconn_2)
    } else {
      if(length(distance_thresholds) > 1){
        names(result_protconn_2) <- paste0("d", distance_thresholds)
      }
    }
    return(result_protconn_2)
  }), error = function(err)err)

  ## End
  if (inherits(result_protconn, "error")){
    stop(result_protconn)
  } else {
    if (length(transboundary) > 1){
      names(result_protconn) <- paste0("Transboundary_", transboundary)
    } else {
      if(lengths(result_protconn) == 1){
        if(isFALSE(plot) & isFALSE(dPC)){
          result_protconn <- result_protconn[[1]][[1]][[1]]
        } else {
          result_protconn <- result_protconn[[1]][[1]]
        }
      } else {
        result_protconn <- result_protconn[[1]]
      }
    }
  }
  return(result_protconn)
}

