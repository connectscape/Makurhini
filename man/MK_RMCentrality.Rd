% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MK_RMCentrality.R
\name{MK_RMCentrality}
\alias{MK_RMCentrality}
\title{Estimate radial and medial centrality measures.}
\usage{
MK_RMCentrality(
  nodes,
  distance = list(type = "centroid"),
  distance_thresholds = NULL,
  binary = TRUE,
  probability = NULL,
  rasterparallel = FALSE,
  write = NULL,
  intern = TRUE
)
}
\arguments{
\item{nodes}{Object containing nodes (e.g., habitat patches or fragments) information. It can be of the following classes:\cr
-   \code{Data.frame} with at least two columns: the first for node IDs and the second for attributes. If the `restoration` argument is used, the data frame must include a third column for restoration values.\cr
-   Spatial data of type vector (class \code{sf, SpatVector, SpatialPolygonsDataFrame}). It must be in a projected coordinate system.\cr
-   Raster (class \code{RasterLayer, SpatRaster}). It must be in a projected coordinate system. The values must be integers representing the ID of each habitat patch or node, with non-habitat areas represented by NA values (see \link[raster]{clump} or \link[terra]{patches}).}

\item{distance}{A \code{matrix} or \code{list} to establish the distance between each pair of nodes. Distance between nodes may be Euclidean distances (straight-line distance) or effective distances (cost distances) by considering the landscape resistance to the species movements. If it is a matrix, then the number of columns and rows must be equal to the number of nodes. This distance matrix could be generated by the \link[Makurhini]{distancefile} function.  If it is a list, then it must contain the distance parameters necessary to calculate the distance between nodes. For example, two of the most important parameters: \code{“type”} and \code{“resistance”}. For \code{"type"} choose one  of the distances:  \bold{"centroid" (faster), "edge", "least-cost" or "commute-time"}. If the type is equal to \code{"least-cost"} or \code{"commute-time"}, then you must use the \code{"resistance"} argument. To see more arguments see the \link[Makurhini]{distancefile} function.}

\item{distance_thresholds}{A \code{numeric} indicating the dispersal distance or distances (meters) of the considered species. If \code{NULL} then distance is estimated as the median dispersal distance between nodes. Alternatively, the \link[Makurhini]{dispersal_distance} function can be used to estimate the dispersal distance using the species home range.}

\item{binary}{\code{logical}. Binary metrics, it only considers the distance thresholds to establish if a pair of nodes is (1) or not connected (0). Probability argument is not necessary.}

\item{probability}{A \code{numeric} value indicating the probability that corresponds to the distance specified in the \code{distance_threshold}. For example, if the \code{distance_threshold} is a median dispersal distance, use a probability of 0.5 (50\%). If the \code{distance_threshold} is a maximum dispersal distance, set a probability of 0.05 (5\%) or 0.01 (1\%). If \code{probability = NULL}, then a probability of 0.5 will be used.}

\item{rasterparallel}{\code{logical}. If nodes is "raster" then you can use this argument to assign the metrics values to the nodes raster. It is useful when raster resolution is less than 100 m<sup>2</sup>.}

\item{write}{\code{character}. Write output \code{sf} object. It is necessary to specify the "Folder direction"
+ "Initial prefix",  for example, \code{"C:/ejemplo"}.}

\item{intern}{\code{logical}. Show the progress of the process, \code{default = TRUE}. Sometimes the advance process does not reach 100 percent when operations are carried out very quickly.}
}
\description{
Use this function to calculate radial and medial centrality measure under one or several distance thresholds.
}
\details{
This function implements Patch-Scale Connectivity or Centrality Measures. Radial measures: degree, strength (using probability argument, for weighted graphs),
eigenvector centrality (eigen), and closeness centrality (close). Medial measures: betweenness centrality (BWC),
node memberships (cluster), and modularity (modules, using probability argument).
The function builds on functions out of Csardi’s ’igraph’ package.
}
\examples{
\dontrun{
library(Makurhini)
library(sf)
data("habitat_nodes", package = "Makurhini")
nrow(habitat_nodes) # Number of patches
#Two distance threshold,
centrality_test <- MK_RMCentrality(nodes = habitat_nodes,
                                distance = list(type = "centroid"),
                                 distance_thresholds = c(10000, 100000),
                                 probability = 0.05,
                                 write = NULL)
centrality_test
plot(centrality_test$d10000["degree"], breaks = "jenks")
plot(centrality_test$d10000["BWC"], breaks = "jenks")
plot(centrality_test$d10000["cluster"], breaks = "jenks")
plot(centrality_test$d10000["modules"], breaks = "jenks")
}
}
\references{
Borgatti, S. P., & Everett, M. G. (2006). A Graph-theoretic perspective on centrality. Social Networks, 28(4), 466–484. https://doi.org/10.1016/j.socnet.2005.11.005\cr
Hanski, I. and Ovaskainen, O. 2000. The metapopulation capacity of a fragmented landscape. Nature 404: 755–758.
}
