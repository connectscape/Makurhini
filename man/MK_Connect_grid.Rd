% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MK_Connect_grid.R
\name{MK_Connect_grid}
\alias{MK_Connect_grid}
\title{Connectivity indexes in a regular grid}
\usage{
MK_Connect_grid(
  nodes,
  area_unit = "ha",
  region = NULL,
  grid = list(hexagonal = TRUE, cellsize = NULL, grid_boundary = FALSE, clip = FALSE,
    tolerance = NULL),
  protconn = TRUE,
  distance_threshold = NULL,
  probability = NULL,
  transboundary = NULL,
  distance = list(type = "centroid"),
  intern = TRUE,
  parallel = NULL
)
}
\arguments{
\item{nodes}{\code{sf, SpatVector, SpatialPolygonsDataFrame}. Object containing nodes (e.g., habitat patches or fragments) of each time to analyze information. Nodes are spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.}

\item{area_unit}{\code{character}. (\emph{optional, default = } \code{"m2"}) \cr. A \code{character} indicating the area units when \code{attribute} is \code{NULL}. Some options are "m2" (the default), "km2", "cm2", or "ha";  See \link[Makurhini]{unit_convert} for details.}

\item{region}{object of class \code{sf}, \code{SpatialPolygonsDataFrame}. Polygon delimiting the region or study area. It must be
in a projected coordinate system.}

\item{grid}{\code{list} or object of class \code{sf}, \code{SpatialPolygonsDataFrame}.
Use this parameter to generate a grid indicating its characteristics in a \code{list} (see \link[Makurhini]{get_grid}) or enter the name of an sf class \code{sf} or \code{SpatialPolygonsDataFrame} with the grid whose coordinate system must be the same as that of the \code{nodes}.
 Example for generating 100 km\out{<sup>2</sup>} hexagons:\cr
\code{list(hexagonal = TRUE, cellsize = unit_convert(100, "km2", "m2"), grid_boundary = FALSE, clip = FALSE, tolerance = NULL)}.}

\item{protconn}{\code{logical}. If \code{TRUE} then the \bold{ProtConn} will be estimated; otherwise, the \bold{PC} index will be estimated.}

\item{distance_threshold}{A \code{numeric} indicating the dispersal distance (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range.}

\item{probability}{A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). Use in case of selecting the \code{"PC"} metric. If \code{probability = NULL}, then a probability of 0.5 will be used.}

\item{transboundary}{\code{numeric}. Buffer to select transboundary polygones, see \link[Makurhini]{MK_ProtConn}.}

\item{distance}{A \code{list} of parameters to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. \cr
 This list must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. For example: \code{distance(type = "least-cost", resistance = raster_resistance)}. \cr
To see more arguments see the \link[Makurhini]{distancefile} function.}

\item{intern}{\code{logical}. Show the progress of the process, \code{default = TRUE}. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.}

\item{parallel}{\code{numeric}. Specify the number of cores to use for parallel processing, \code{default = NULL}. Parallelize the function using furrr package.}
}
\description{
Use the function to compute the Protected Connected (ProtConn), EC, PC or IIC indexes in a regular grid.
}
\examples{
\dontrun{
library(Makurhini)
library(sf)
load(system.file("extdata", "Protected_areas.rda",
                package = "Makurhini", mustWork = TRUE))
data("Ecoregions", package = "Makurhini")
ecoregion <- Ecoregions[1,]
plot(st_geometry(ecoregion), col = "#7E6A9F")
#ProtConn
hexagons_priority <- MK_Connect_grid(nodes = Protected_areas,
                                    region = ecoregion,
                                    area_unit = "ha",
                                    grid = list(hexagonal = TRUE,
                                                cellsize = unit_convert(5000, "km2", "m2")),
                                    protconn = TRUE,
                                    distance_threshold = 3000,
                                    probability = 0.5,
                                    transboundary = 6000,
                                    distance = list(type = "centroid"),
                                    intern = TRUE,
                                    parallel = NULL)
hexagons_priority
plot(st_geometry(ecoregion), col = "#7E6A9F")
plot(hexagons_priority["ProtConn"], add = TRUE)
}
}
\references{
Matt Strimas-Mackey. \url{http://strimas.com/spatial/hexagonal-grids/}.\cr
Saura, S., Bastin, L., Battistella, L., Mandrici, A., & Dubois, G. (2017). Protected areas in the
world's ecoregions: How well connected are they? Ecological Indicators, 76, 144–158.
Saura, S. & Torne, J. (2012). Conefor 2.6. Universidad Politecnica de Madrid. Available
 at \url{www.conefor.org}.\cr
Pascual-Hortal, L. & Saura, S. (2006). Comparison and development of new graph-based landscape
 connectivity indices: towards the priorization of habitat patches and corridors for conservation.
 Landscape Ecology, 21(7): 959-967.\cr
Saura, S. & Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity
in landscape conservation planning: comparison with existing indices and application to a case study.
Landscape and Urban Planning, 83(2-3): 91-103.
}
