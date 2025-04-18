% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MK_Fragmentation.R
\name{MK_Fragmentation}
\alias{MK_Fragmentation}
\title{Fragmentation Statistics}
\usage{
MK_Fragmentation(
  nodes = NULL,
  edge_distance = 500,
  min_node_area = 100,
  landscape_area = NULL,
  area_unit = "ha",
  perimeter_unit = "km",
  plot = FALSE,
  write = NULL
)
}
\arguments{
\item{nodes}{Object of class \code{sf, sfc, sfg, SpatialPolygons}. Individual nodes, the object must be in a projected coordinate system.}

\item{edge_distance}{\code{numeric}. Distance to edge in meters. Default equal 500 m (Haddad et al. 2015)}

\item{min_node_area}{\code{numeric}. Minimum node area used to calculate the number of nodes with an area smaller than the one provided. Default equal 100 km\out{<sup>2</sup>} (Haddad et al. 2015)}

\item{landscape_area}{\code{numeric}. Total area of the study landscape in km\out{<sup>2</sup>} (optional). If NULL the total nod area will be used.}

\item{area_unit}{\code{character}. You can set an area unit (e.g., "km2", "cm2", "m2", "ha"; see \link[Makurhini]{unit_convert}). Default equal to square kilometers "km2".}

\item{perimeter_unit}{\code{character}. You can set a perimeter unit (e.g., "km", "cm", "m", "ha"; see \link[Makurhini]{unit_convert}). Default equal to kilometers "km".}

\item{plot}{\code{logical}. Basic histograms and core area - edge map.}

\item{write}{\code{character}. Write the table (landscape statistics), sf object (patch/node statistics) and plots. It's necessary to specify the path and prefix, for example,
to save in the path "C:/Folder" with the prefix "Fragmentation": \code{"C:/Folder/Fragmentation"}.}
}
\value{
Patch/node and landscape statistics:\cr
1) Patches Area in square kilometers.\cr
2) Number patches.\cr
3) Mean size of patches.\cr
4) Number of patches smaller than the parameter \code{min_node_area} (default = 100 km\out{<sup>2</sup>}).\cr
5) Percentage of patches smaller than the parameter \code{min_node_area} (default = 100 km\out{<sup>2</sup>}).\cr
6) Total patch edge. Total perimeter of the patches (unit = \code{perimeter_unit}).\cr
7) Edge density. Total perimeter per unit of area (unit = \code{area_unit}), default = km\out{<sup>2</sup>}. A value of 0 is present when there is no edge in the landscape.\cr
8) Patch density.\cr
9) Total core area (units = \code{area_unit}) considering the distance set in the parameter \code{edge_distance} (delfault = 500 m).\cr
10) Cority index. It is a measure of fragmentation with respect to a distance from the core area (parameter \code{edge_distance}; delfault = 500 m), where a value of 1 indicates a landscape without fragmentation. Average for landscape level.\cr
11) Shape Index. A simple shape metric that takes values from 1 (perfectly compact) to infinity is derived by dividing the perimeter by the perimeter of a circle of the same area. Average for landscape level. \cr
12) Fractal dimension. The index reflects the complexity of the shape of the fragment. A fractal dimension greater than 1 indicates an increase in the complexity of the shape. When the value is close to 1 the shape is simple, such as squares.Average for landscape level. \cr
13) Effective Mesh Size. Effective Mesh Size (MESH) is a measure of the degree of fragmentation in the landscape ranging from 0 to the total landscape area. MESH is maximum when the landscape unit consists of a single habitat fragment or the habitat is continuous beyond the landscape unit analyzed (Moser, 2007).\cr
14) Core percent (patch level). Percentage of core area in the patch (units = \code{area_unit}) considering the distance set in the parameter \code{edge_distance} (delfault = 500 m).\cr
15) Edge percent (patch level). Percentage of edge in the patch (units = \code{area_unit}) considering the distance set in the parameter \code{edge_distance} (delfault = 500 m).\cr
16) PARA (patch level). Ratio of the patch perimeter to area.
*NOTE.* In the results we use the term patches instead of nodes due to the common use of this term in fragmentation statistics in science.
}
\description{
Calculate patches/nodes and landscape statistics
}
\examples{
data("habitat_nodes", package = "Makurhini")
nrow(habitat_nodes) # Number of nodes
fragmentation <- MK_Fragmentation(nodes = habitat_nodes, edge_distance = 1000, plot = TRUE)
#Table
fragmentation$`Summary landscape metrics (Viewer Panel)`
#Shapefile
fragmentation$`Patch statistics shapefile`
}
\references{
Haddad et al. (2015). Science Advances 1(2):e1500052. https://www.science.org/doi/10.1126/sciadv.1500052.\cr
McGarigal, K., S. A. Cushman, M. C. Neel, and E. Ene. 2002. FRAGSTATS: Spatial Pattern Analysis Program for Categorical Maps. Computer software program produced by the authors at the University of Massachusetts, Amherst. Available at the following web site:
 \url{www.umass.edu/landeco/research/fragstats/fragstats.html}.\cr
Moser, B., Jaeger, J.A.G., Tappeiner, U. et al. Modification of the effective mesh size for measuring landscape fragmentation to solve the boundary problem. Landscape Ecol 22, 447–459 (2007).  \url{https://doi.org/10.1007/s10980-006-9023-0}
}
