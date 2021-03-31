#' @title Cores shapefile for analysis
#'
#' @description A simple feature collection \code{sf} with 142 Vegetation patches with less alteration or degradation that
#'  is for testing and running of MK_fragmentation and MK_dPCIIC examples.
#' @format A \code{sf} with 142 features and 1 field, which are:
#' \describe{
#' \item{id}{Vegetation patches identifier.}
#'  }
#' @references
#' INEGI. (2013). Conjunto de datos vectoriales de uso del suelo y vegetacion, serie V (capa union), escala 1:250, 000.
#' Instituto Nacional de Estadistica y Geografia, Aguascalientes.
"vegetation_patches"


#' @title List of vegetation patches
#' @description A list with 4 shapefiles. Each shapefile correspond to vegetation patches with
#' less alteration or degradation  in the region of the Biosphere Reserve Mariposa Monarca for
#' the nex periods: "1993", "2003", "2007", "2011". The list is used for testing and running
#' of MK_dECA and test_ECA_distance examples.
#' @format A \code{sf} with 4 shapefiles:
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
#' @format A \code{sf} with 4 shapefiles:
"study_area"
