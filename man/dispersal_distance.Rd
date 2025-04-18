% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dispersal_distance.R
\name{dispersal_distance}
\alias{dispersal_distance}
\title{Estimate dispersal distance using home range}
\usage{
dispersal_distance(
  home_range = NULL,
  home_range_unit = "ha",
  dispersal_distance_unit = "m"
)
}
\arguments{
\item{home_range}{\code{numeric}. Species home range}

\item{home_range_unit}{\code{character}. You can set an area unit (e.g., "km2", "cm2", "ha"; see \link[Makurhini]{unit_convert}).
Default equal to hectares \code{"ha"}.}

\item{dispersal_distance_unit}{\code{character}. You can set an distance unit (e.g., "km", "m", "cm"; see \link[Makurhini]{unit_convert}).
Default equal to meters \code{"m"}.}
}
\value{
Returns a class object \code{data.frame} containing the scattering distances for a given home range.
}
\description{
Estimate the dispersal distance using the home range and the equation
proposed by Bowman, Jaeger and Fahrig (2002)
}
\examples{
dispersal_distance(home_range = 0.6,
                   home_range_unit = "ha",
                   dispersal_distance_unit = "m")
}
\references{
Bowman, J., Jaeger, J. A., & Fahrig, L. (2002). Dispersal distance of mammals is proportional to home range size. Ecology, 83(7), 2049-2055.
}
