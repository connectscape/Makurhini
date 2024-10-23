## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(Makurhini)
library(sf)
library(raster)
library(ggplot2)
library(classInt)
library(ggspatial)
library(viridis)

## ----message=FALSE, warning=FALSE---------------------------------------------
#Habitat nodes
data("habitat_nodes", package = "Makurhini")
nrow(habitat_nodes)

#Study area
data("TMVS", package = "Makurhini")

#Resistance
data("resistance_matrix", package = "Makurhini")


## ----message=FALSE, warning=FALSE---------------------------------------------
raster_map <- as(resistance_matrix, "SpatialPixelsDataFrame")
raster_map <- as.data.frame(raster_map)
colnames(raster_map) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data = raster_map, aes(x = x, y = y, fill = value), alpha = 0.8) + 
  geom_sf(data = TMVS, aes(color = "Study area"), fill = NA, color = "black") +
  geom_sf(data = habitat_nodes, aes(color = "Habitat nodes"), fill = "forestgreen") +
  scale_fill_viridis(option = "B", name = "Resistance")+
  scale_color_manual(name = "", values = "black")+
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  annotation_scale(
    location = "bl",
    bar_cols = c("grey10", "white"),
    text_family = "ArcherPro Book"
  ) +
  annotation_north_arrow(
    location = "br", which_north = "true",
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_orienteering(
      fill = c("grey10", "white"),
      line_col = "grey5",
      text_family = "ArcherPro Book"
    )
  )

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  PC_example_1 <- MK_dPCIIC(nodes = habitat_nodes,
#                          attribute = NULL,
#                          distance = list(type = "centroid"),
#                          parallel = NULL,
#                          metric = "PC",
#                          probability = 0.5,
#                          distance_thresholds = c(250, 1500, 3000, 10000))

## ----eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE----------------------
PC_example_1 <- readRDS("E:/Makurhini/Paper_Makurhini/Ejemplos/PC_example_1.rds")

## ----eval=TRUE, message=FALSE, warning=FALSE----------------------------------
class(PC_example_1)

names(PC_example_1)

head(PC_example_1$d10000)

## ----eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE----------------------
interv <- classIntervals(PC_example_1$d10000$dPC, 9, "jenks")[[2]] #9 intervalos
ggplot()+
  geom_sf(data = TMVS)+
  geom_sf(data = PC_example_1$d10000, aes(fill = cut(dPC, breaks = interv)), color = NA)+
  scale_fill_brewer(type = "qual",
                    palette = "RdYlGn",
                    name = "dPC",
                    na.translate = FALSE)+
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1,0.21),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.text = element_text(size = 5.5, family = "Times"),
    legend.title = element_text(size = 5.5, family = "Times")
  ) + labs(title="Euclidean distance")


## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  PC_example_2 <- MK_dPCIIC(nodes = habitat_nodes,
#                          attribute = NULL,
#                          distance = list(type = "least-cost",
#                                          resistance = resistance_matrix),
#                          parallel = NULL,
#                          metric = "PC",
#                          probability = 0.5,
#                          distance_thresholds = c(250, 1500, 3000, 10000))

## ----eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE----------------------
PC_example_2 <- readRDS("E:/Makurhini/Paper_Makurhini/Ejemplos/PC_example_2.rds")

## ----eval=TRUE, message=FALSE, warning=FALSE----------------------------------
class(PC_example_2)

names(PC_example_2)

head(PC_example_2$d10000)

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  write_sf(PC_example$d10000, “.../dPC_d0000.shp”)

## ----eval=TRUE, message=FALSE, warning=FALSE----------------------------------
#Keep the same range of values of PC_example_1 for comparison, only the highest range changes.
interv[length(interv)] <- max(PC_example_2$d10000$dPC)
ggplot()+
  geom_sf(data = TMVS)+
  geom_sf(data = PC_example_2$d10000, aes(fill = cut(dPC, breaks = interv)), color = NA)+
  scale_fill_brewer(type = "qual",
                    palette = "RdYlGn",
                    name = "dPC",
                    na.translate = FALSE)+
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.21),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.text = element_text(size = 5.5, family = "Times"),
    legend.title = element_text(size = 5.5, family = "Times")
  )+ labs(title="Least-cost distance")


## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  PC_example <- MK_dPCIIC(nodes = patches, attribute = NULL, area_unit = "ha",
#                          distance = list(type = "least-cost",
#                                          resistance = resistance_matrix,
#                                          least_cost.java = TRUE,
#                                          cores.java = 4),
#                          metric = "PC", overall = TRUE, probability = 0.5,
#                          distance_thresholds = c(250, 1500, 3000, 10000),
#                          intern = FALSE)
#  #155.59 sec elapsed

## ----eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE----------------------
PC_example_3 <- MK_dPCIIC(nodes = habitat_nodes,
                        attribute = NULL,
                        area_unit = "ha",
                        distance = list(type = "centroid"),
                        parallel = NULL,
                        metric = "PC",
                        probability = 0.5,
                        distance_thresholds = c(250, 1500, 3000, 10000),
                        overall = TRUE, intern = FALSE)

## ----eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE----------------------
#  PC_example_3 <- MK_dPCIIC(nodes = habitat_nodes,
#                          attribute = NULL,
#                          area_unit = "ha",
#                          distance = list(type = "centroid"),
#                          parallel = NULL,
#                          metric = "PC",
#                          probability = 0.5,
#                          distance_thresholds = c(250, 1500, 3000, 10000),
#                          overall = TRUE)
#  

## -----------------------------------------------------------------------------
class(PC_example_3)
class(PC_example_3$d10000$node_importances_d10000)
class(PC_example_3$d10000$overall_d10000)
PC_example_3$d10000$overall_d10000

## -----------------------------------------------------------------------------
PC_example_3$d10000$node_importances_d10000

## -----------------------------------------------------------------------------
PC_example_3$d10000$overall_d10000

## ----eval=TRUE, message=FALSE, warning=FALSE, echo=FALSE----------------------
#Maximum landcape attribute or LA = total area of the estudy area
Area <- unit_convert( st_area(TMVS), "m2", "ha") #hectares
PC_example_3 <- MK_dPCIIC(nodes = habitat_nodes,
                        attribute = NULL,
                        area_unit = "ha",
                        distance = list(type = "centroid"),
                        parallel = NULL,
                        metric = "PC",
                        probability = 0.5,
                        distance_thresholds = c(250, 1500, 3000, 10000),
                        LA = Area,
                        onlyoverall = TRUE, intern = FALSE)

## ----eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE----------------------
#  #Maximum landcape attribute or LA = total area of the estudy area
#  Area <- unit_convert( st_area(TMVS), "m2", "ha") #hectares
#  PC_example_3 <- MK_dPCIIC(nodes = habitat_nodes,
#                          attribute = NULL,
#                          area_unit = "ha",
#                          distance = list(type = "centroid"),
#                          parallel = NULL,
#                          metric = "PC",
#                          probability = 0.5,
#                          distance_thresholds = c(250, 1500, 3000, 10000),
#                          LA = Area,
#                          onlyoverall = TRUE)

## ----eval=TRUE----------------------------------------------------------------
class(PC_example_3)
class(PC_example_3$d10000)
PC_example_3$d10000

## ----eval=TRUE, echo=FALSE, message=FALSE, warning=FALSE----------------------
IIC_example_4 <- MK_dPCIIC(nodes = habitat_nodes,
                        attribute = NULL,
                        area_unit = "ha",
                        distance = list(type = "centroid"),
                        parallel = NULL,
                        metric = "IIC",
                        probability = NULL,
                        distance_thresholds = c(250, 1500, 3000, 10000),
                        LA = Area,
                        overall = TRUE, intern = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  IIC_example_4 <- MK_dPCIIC(nodes = habitat_nodes,
#                          attribute = NULL,
#                          area_unit = "ha",
#                          distance = list(type = "centroid"),
#                          parallel = NULL,
#                          metric = "IIC",
#                          probability = NULL,
#                          distance_thresholds = c(250, 1500, 3000, 10000),
#                          LA = Area,
#                          overall = TRUE)

## ----eval=TRUE----------------------------------------------------------------
class(IIC_example_4)
class(IIC_example_4$d10000$node_importances_d10000)
class(IIC_example_4$d10000$overall_d10000)

## -----------------------------------------------------------------------------
IIC_example_4$d10000$node_importances_d10000

## -----------------------------------------------------------------------------
IIC_example_4$d10000$overall_d10000

