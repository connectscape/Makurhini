## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", message=FALSE, warning=FALSE
)

## ----setup--------------------------------------------------------------------
library(Makurhini)
library(sf)
library(ggplot2)

## ----polygons-----------------------------------------------------------------
data("vegetation_patches", package = "Makurhini")
nrow(vegetation_patches) # Number of patches
plot(st_geometry(vegetation_patches), col = "#00B050")

## -----------------------------------------------------------------------------
Fragmentation_test <- MK_Fragmentation(patches = vegetation_patches, edge_distance = 500, plot = F, min_patch_area = 100, landscape_area = NULL, area_unit = "km2", perimeter_unit = "km")

## -----------------------------------------------------------------------------
names(Fragmentation_test)
Fragmentation_test$`Summary landscape metrics (Viewer Panel)`

## -----------------------------------------------------------------------------
head(Fragmentation_test[[2]])

## ----echo=FALSE---------------------------------------------------------------
plot(Fragmentation_test[[2]]["CAPercent"])

## ----echo=FALSE---------------------------------------------------------------
plot(Fragmentation_test[[2]]["EdgePercent"])

## ----echo=FALSE---------------------------------------------------------------
plot(Fragmentation_test[[2]]["PARA"])

## ----echo=FALSE---------------------------------------------------------------
plot(Fragmentation_test[[2]]["ShapeIndex"])

## ----echo=FALSE---------------------------------------------------------------
plot(Fragmentation_test[[2]]["FRAC"])

## ----echo=FALSE---------------------------------------------------------------
library(purrr)
Fragmentation_test.2 <- map_dfr(seq(100, 1000, 100), function(x){
  x.1 <- MK_Fragmentation(patches = vegetation_patches, edge_distance = x, plot = FALSE)[[2]]
  CA <- mean(x.1$CAPercent)
  Edge <- mean(x.1$EdgePercent)
  x.2 <- rbind(data.frame('Edge distance' = x, Type = "Core Area", Percentage = CA),
                     data.frame('Edge distance' = x, Type = "Edge", Percentage = Edge))
  return(x.2)
})

head(Fragmentation_test.2)

## ----echo=FALSE---------------------------------------------------------------
ggplot(Fragmentation_test.2, aes(x = Edge.distance, y = Percentage, group = Type)) +
  geom_line(aes(color = Type))+
  geom_point(aes(color = Type))+ ylim(0,100)+
  scale_x_continuous("Edge depth distance (m)", labels = as.character(Fragmentation_test.2$Edge.distance), breaks = Fragmentation_test.2$Edge.distance)+
  scale_color_brewer(palette="Dark2")+
  theme_classic()

