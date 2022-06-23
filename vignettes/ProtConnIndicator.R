## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", message=FALSE, warning=FALSE
)

## ----setup--------------------------------------------------------------------
library(Makurhini)
library(raster)
library(sf)

## ----polygons-----------------------------------------------------------------
data("Protected_areas", package = "Makurhini")
data("regions", package = "Makurhini")

## ----echo=FALSE---------------------------------------------------------------
plot(regions, col=c("#FC8D62", "#8DA0CB", "#E78AC3"))
plot(Protected_areas, col = "#00B050", add=TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
region <- regions[1,]
test.1 <- MK_ProtConn(nodes = Protected_areas, region = region, 
                      area_unit = "ha", distance = list(type= "centroid"), 
                      distance_thresholds = c(10000, 30000), probability = 0.5, 
                      transboundary = 50000, transboundary_type = "region",
                      LA = NULL, plot = TRUE, write = NULL,
                      intern = FALSE)

## -----------------------------------------------------------------------------
test.1$d10000$`Protected Connected (Viewer Panel)`
test.1$d10000

## -----------------------------------------------------------------------------
test.1$d30000$`Protected Connected (Viewer Panel)`
test.1$d30000

## ----message=FALSE, warning=FALSE---------------------------------------------
test.2 <- MK_ProtConnMult(nodes = Protected_areas, regions = regions, 
                          area_unit = "ha", distance = list(type= "centroid"), 
                          distance_thresholds = c(10000, 30000), probability = 0.5, 
                          transboundary = 50000, transboundary_type = "region",
                          plot = FALSE, write = NULL, 
                          parallel = NULL, intern = FALSE)

## -----------------------------------------------------------------------------
names(test.2)
test.2$ProtConn_10000$ProtConn_overall10000

## -----------------------------------------------------------------------------
test.2$ProtConn_10000$ProtConn_10000

## ----echo=FALSE---------------------------------------------------------------
ProtConn <- test.2[[1]]$ProtConn_10000
plot(ProtConn["ProtConn"])

