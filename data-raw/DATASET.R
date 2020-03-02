## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")

ruta <- system.file("extdata", "Fragmentation.RData", package = "Makurhini")
load(ruta)
use_data(vegetation_patches, compress = "xz")

ruta <- system.file("extdata", "WDPA_May2019_MEX-shapefile-polygons.shp", package = "Makurhini")
Protected_areas <- shapefile(ruta)
Protected_areas <- Protected_areas[,1]

ruta <- system.file("extdata", "Ecoregions2017.shp", package = "Makurhini")
regions <- shapefile(ruta)
regions <- regions[,1]
use_data(Protected_areas, regions, compress = "xz", overwrite = T)

use_data(vegetation_patches, raster_vegetation_patches, compress = "xz", overwrite = T)
