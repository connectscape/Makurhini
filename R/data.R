#' @title Nodes shapefile for analysis
#'
#' @description A simple feature collection \code{sf} with 404 habitat patches/nodes with less alteration or degradation that
#'  is for testing and running of MK_fragmentation and MK_dPCIIC examples.
#' @format A \code{sf} with 404 features and 1 field.
#' \describe{
#' \item{id}{Vegetation patches identifier.}
#'  }
#' @references
#' Correa Ayram, C. A., Mendoza, M. E., Etter, A., & Pérez Salicrup, D. R. (2017). Anthropogenic impact on habitat connectivity: A multidimensional human footprint index evaluated in a highly biodiverse landscape of Mexico. Ecological Indicators, 72, 895-909. https://doi.org/10.1016/j.ecolind.2016.09.007
"habitat_nodes"


#' @title List of vegetation patches
#' @description A list with 4 shapefiles. Each shapefile correspond to vegetation patches with
#' less alteration or degradation  in the region of the Biosphere Reserve Mariposa Monarca for
#' the nex periods: "1993", "2003", "2007", "2011". The list is used for testing and running
#' of MK_dECA and test_ECA_distance examples.
#' @format A \code{sf} with 4 shapefiles.
#' @references
#' INEGI. (2001). Conjunto de datos vectoriales de uso del suelo y vegetacion, serie II (Continuo nacional), escala 1:250, 000.
#' Instituto Nacional de Estadistica y Geografia, Aguascalientes.
#' INEGI. (2005). Conjunto de datos vectoriales de uso del suelo y vegetacion, serie III (Continuo nacional), escala 1:250, 000.
#' Instituto Nacional de Estadistica y Geografia, Aguascalientes.
#' INEGI. (2009). Conjunto de datos vectoriales de uso del suelo y vegetacion, serie IV (Continuo nacional), escala 1:250, 000.
#' Instituto Nacional de Estadistica y Geografia, Aguascalientes.
#' INEGI. (2013). Conjunto de datos vectoriales de uso del suelo y vegetacion, serie V (capa union), escala 1:250, 000.
#' Instituto Nacional de Estadistica y Geografia, Aguascalientes.
"list_forest_patches"


#' @title Study area example
#' @description shapefile. Biosphere reserve Mariposa Monarca plus a buffer of 10 km. It is
#' used to run the MK_dECA and test_ECA_distance examples.
#' @format A \code{sf}
"study_area"


#' @title Trans-Mexican Volcanic System (TMVS)
#' @description shapefile. Trans-Mexican Volcanic System (TMVS).
#' @format A \code{sf}
"TMVS"


#' @title habitat nodes in raster format
#'
#' @description A simple feature collection \code{sf} with 404 habitat patches/nodes with less alteration or degradation that
#'  is for testing and running of MK_fragmentation and MK_dPCIIC examples.
#' @format A \code{raster} with 404 integer values.
#' \describe{
#' \item{id}{habitat patches identifier.}
#'  }
#' @references
#' Correa Ayram, C. A., Mendoza, M. E., Etter, A., & Pérez Salicrup, D. R. (2017). Anthropogenic impact on habitat connectivity: A multidimensional human footprint index evaluated in a highly biodiverse landscape of Mexico. Ecological Indicators, 72, 895-909. https://doi.org/10.1016/j.ecolind.2016.09.007
"habitat_nodes_raster"


#' @title Landscape resistance to dispersal
#'
#' @description A RasterLayer with Human footprint values. The raster was aggregated by a factor of 5 to change its original resolution from 100m to 500m.
#' @format A \code{raster} with Human footprint values from 0 to 100:
#' \describe{
#' \item{values}{Human footprint values from 0 to 100}
#'  }
#' @references
#' Correa Ayram, C. A., Mendoza, M. E., Etter, A., & Pérez Salicrup, D. R. (2017). Anthropogenic impact on habitat connectivity: A multidimensional human footprint index evaluated in a highly biodiverse landscape of Mexico. Ecological Indicators, 72, 895-909. https://doi.org/10.1016/j.ecolind.2016.09.007
"resistance_matrix"


#' @title Protected areas raster for analysis
#'
#' @description raster \code{raster} with 706 polygons of Protected Areas in Mexico.
#' @references
#' UNEP-WCMC (2019). Protected Area Profile for Mexico from the World Database of Protected Areas, May 2019.
#' Available at: www.protectedplanet.net
"Protected_areas_raster"


#' @title Ecoregions shapefile for analysis
#'
#' @description A Spatial Polygons DataFrame \code{SpatialPolygonsDataFrame} with 3 ecoregions of Mexico.
#' @format A \code{SpatialPolygonsDataFrame} with 3 features and 3 fields, which are:
#' \describe{
#' \item{OBJECTID}{Ecoregions identifier.}
#'  }
#' @references
#'Dinerstein, E., Olson, D., Joshi, A., Vynne, C., Burgess, N. D., Wikramanayake, E., and Kindt, R. (2017).
#'An Ecoregion-Based Approach to Protecting Half the Terrestrial Realm. BioScience, 67(6), 534?545. https://doi.org/10.1093/biosci/bix014
#'https://ecoregions2017.appspot.com
"Ecoregions"

#' @title File of habitat nodes in Spain
#'
#' @description A \code{data.frame} with 359 habitat patches/nodes identified as
#' the areas covered by forest according to the CORINE land cover map from Spain.
#' @format A \code{data.frame} with 359 files and 2 columns.
#' \describe{
#' \item{ID}{nodes identifier}
#' \item{attribute}{area in ha}
#' }
#' @references
#' European Commission. (2020). COMMUNICATION FROM THE COMMISSION TO THE EUROPEAN PARLIAMENT, THE COUNCIL, THE EUROPEAN ECONOMIC AND SOCIAL COMMITTEE AND THE COMMITTEE OF THE REGIONS EU Biodiversity Strategy for 2030 Bringing nature back into our lives. https://ec.europa.eu/research/environment/index.cfm?pg=nbs \cr
#' Goicolea, T., Cisneros-Araújo, P., Vega, C.A. et al. Landscape connectivity for predicting the spread of ASF in the European wild boar population. Sci Rep 14, 3414 (2024). https://doi.org/10.1038/s41598-024-53869-5
"habitat_nodes_spain"


#' @title Distance matrix between habitat nodes in Spain
#'
#' @description A \code{Large matrix} with 128881 elements or 359 files and 359 columns. The distance file gives the current effective distance of corridors between patches.
#' @format A \code{matrix} with 359 files and 359 columns.
#' @references
#' Goicolea, T., Cisneros-Araújo, P., Vega, C.A. et al. Landscape connectivity for predicting the spread of ASF in the European wild boar population. Sci Rep 14, 3414 (2024). https://doi.org/10.1038/s41598-024-53869-5
"dist_original"


#' @title Distance matrix between habitat nodes in Spain after an hypothetical restoration
#'
#' @description A \code{Large matrix} with 128881 elements or 359 files and 359 columns. The distance file gives the potential (after hypothetical restoration measures) effective distance of corridors between patches.
#' @format A \code{matrix} with 359 files and 359 columns.
#' @references
#' Goicolea, T., Cisneros-Araújo, P., Vega, C.A. et al. Landscape connectivity for predicting the spread of ASF in the European wild boar population. Sci Rep 14, 3414 (2024). https://doi.org/10.1038/s41598-024-53869-5
"dist_restoration"
