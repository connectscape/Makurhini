## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(Makurhini)
library(sf)
library(ggplot2)
library(rmapshaper)

## ----message=FALSE, warning=FALSE---------------------------------------------
#Protected areas
load(system.file("extdata", "Protected_areas.rda",
                 package = "Makurhini", mustWork = TRUE))
nrow(Protected_areas)

#Ecoregions
data("Ecoregions", package = "Makurhini")
nrow(Ecoregions)
#For practicality we will use the first three columns.
Ecoregions <- Ecoregions[,1:3]


## ----warning=FALSE, message=FALSE---------------------------------------------
mask_ecoregions <- ms_dissolve(Ecoregions)
PAs_national <- ms_clip(Protected_areas, mask_ecoregions)
PAs_transnational <- ms_erase(Protected_areas, mask_ecoregions)
PAs_transnational$Type <- "PAs in neighboring countries"
PAs_subnational <- PAs_national[PAs_national$ESCALA_2 == "Subnacional",]
PAs_subnational$Type <- "Subnational PAs"
PAs_national <- PAs_national[PAs_national$ESCALA_2 == "Nacional",]
PAs_national$Type <- "National PAs"
PAs <- rbind(PAs_national, PAs_subnational, PAs_transnational)
PAs$Type <- factor(PAs$Type, levels = c("National PAs", "Subnational PAs", "PAs in neighboring countries"))

ggplot() +
  geom_sf(data = Ecoregions, aes(fill = "Ecoregions"), color = "black") +
  geom_sf(data = PAs, aes(fill=Type), color = NA) +
  scale_fill_manual(name = "Type", values = c("#1DAB80", "#FF00C5", "#E06936", "#8D8BBE"))+
  theme_minimal() 

## ----echo=FALSE---------------------------------------------------------------
ProtConn_1 <- readRDS("C:/Users/Usuario/Documents/R/TEST_Folder/ProtConn_1.rds")

## ----eval= FALSE, message=FALSE, warning=FALSE--------------------------------
#  #Select first ecoregion
#  Ecoregion_1 <- Ecoregions[1,]
#  
#  #keep = 0.6 simplify the geometry and reduce the number of vertices
#  ProtConn_1 <- MK_ProtConn(nodes = Protected_areas, region = Ecoregion_1,
#                            area_unit = "ha",
#                            distance = list(type= "edge", keep = 0.6),
#                            distance_thresholds = 10000, probability = 0.5,
#                            transboundary = 50000, plot = TRUE, intern = FALSE)
#  

## -----------------------------------------------------------------------------
class(ProtConn_1)
names(ProtConn_1)

## ----eval= TRUE, message=FALSE, warning=FALSE---------------------------------
ProtConn_1$`Protected Connected (Viewer Panel)`

## ----eval= TRUE, message=FALSE, warning=FALSE---------------------------------
ProtConn_1$`ProtConn Plot`

## ----echo=FALSE---------------------------------------------------------------
ProtConn_1 <- readRDS("C:/Users/Usuario/Documents/R/TEST_Folder/ProtConn_1a.rds")

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  ProtConn_1 <- lapply(1:nrow(Ecoregions), function(x){
#    				x.1 <- MK_ProtConn(nodes = Protected_areas,
#    				                   region = Ecoregions[x,],
#                    area_unit = "ha",
#              			distance = list(type= "edge"),
#              			distance_thresholds = 10000,
#              			probability = 0.5, transboundary = 50000,
#              			LA = NULL, plot = TRUE, intern = FALSE)
#  return(x.1) })
#  class(ProtConn_1)

## ----eval=FALSE, message=FALSE, warning=FALSE---------------------------------
#  Ecoregions$Protconn <- sapply(x.1, function(x){
#      	      	   x.1 <- as.data.frame(x)
#      	      	   x.1 <- x.1[3, 4] #Extract ProtConn value
#      	      	   return(x.1)})
#  head(Ecoregions)

## ----eval=FALSE---------------------------------------------------------------
#  write_sf(Ecoregions, "./Protconn.gpkg")
#  

## ----eval=FALSE---------------------------------------------------------------
#  ProtConn_2 <- MK_ProtConnMult(nodes = Protected_areas,
#                                region = Ecoregions,
#                                area_unit = "ha",
#                                distance = list(type= "edge"),
#                                distance_thresholds = 10000,
#                                probability = 0.5, transboundary = 50000,
#                                plot = TRUE, parallel = 4)

## ----echo=FALSE---------------------------------------------------------------
ProtConn_2 <- readRDS("C:/Users/Usuario/Documents/R/TEST_Folder/ProtConn_1b.rds")

## -----------------------------------------------------------------------------
class(ProtConn_2)
names(ProtConn_2)

## ----eval=TRUE----------------------------------------------------------------
ProtConn_2$ProtConn_10000$ProtConn_overall10000

## ----eval=TRUE----------------------------------------------------------------
ProtConn_2$ProtConn_10000$`ProtConn Plot`

## -----------------------------------------------------------------------------
head(ProtConn_2$ProtConn_10000$ProtConn_10000)

## ----eval=TRUE, message=FALSE, warning=FALSE, echo= FALSE---------------------
interv <- c(0.0701, 1.9375, 4.2690, 6.6786, 10.7244, 17.8158, 25.6303, 41.8570, 45.4735, 97.7425)

## ----eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE----------------------
#  #We can use some package to get intervals for example classInt R Packge:
#  #library(classInt)
#  #interv <- classIntervals(ProtConn_2$ProtConn_10000$ProtConn_10000$ProtConn, 9, "jenks")[[2]]

## ----message=FALSE, warning=FALSE---------------------------------------------
ggplot()+
  geom_sf(data = Ecoregions)+
  geom_sf(data = ProtConn_2$ProtConn_10000$ProtConn_10000, 
          aes(fill = cut(ProtConn, breaks = interv)), color = NA)+
  scale_fill_brewer(type = "qual",
                    palette = "RdYlGn",
                    name = "ProtConn",
                    na.translate = FALSE)+
  theme_minimal() +
  theme(
    legend.position.inside = c(0.1,0.21),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

