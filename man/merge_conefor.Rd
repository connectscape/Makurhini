% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_conefor.R
\name{merge_conefor}
\alias{merge_conefor}
\title{Join a connectivity index (dIIC o dPC, and the "intra",  "dflux", "connector") to core areas}
\usage{
merge_conefor(
  datat = NULL,
  pattern = NULL,
  merge_shape = NULL,
  id = NULL,
  dA = FALSE,
  var = FALSE,
  write = NULL
)
}
\arguments{
\item{datat}{Data table, data frame or direction of node_importances.txt files (character).}

\item{pattern}{character. Prefix (e.g."*AREA.txt" )}

\item{merge_shape}{object of class sf, sfc, sfg or SpatialPolygons. It has to have the same "id" used to estimate the node importance.txt}

\item{id}{character. Column name with the core id}

\item{dA}{logical. If TRUE, the delta attribute and its variance is selected. Default = FALSE.}

\item{var}{logical. If TRUE the metric and fractions variance is reteined.}

\item{write}{character. Write the shapefile, example, "C:/ejemplo/sahapefile.shp".}
}
\description{
Join a connectivity index (dIIC o dPC, and the "intra",  "dflux", "connector") to core areas
}
\keyword{internal}
