#' Estimate dispersal distance using home range
#'
#' Estimate the dispersal distance using the home range and the equation
#' proposed by Bowman, Jaeger and Fahrig (2002)
#' @param home_range numeric. Species home range
#' @param home_range_unit character. You can set an area unit (e.g., "km2", "cm2", "ha"; see Makurhini::unit_convert).
#' Default equal to hectares "ha".
#' @param dispersal_distance_unit character. You can set an distance unit (e.g., "km", "m", "cm"; see Makurhini::unit_convert).
#' Default equal to meters "m".
#' @param intern logical. Show the result, default = TRUE
#' @references Bowman, J., Jaeger, J. A., & Fahrig, L. (2002). Dispersal distance of mammals is proportional to home range size. Ecology, 83(7), 2049-2055.
#' @export
#' @examples
#' dispersal_distance(home_range = 0.6,
#'                    home_range_unit = "ha",
#'                    dispersal_distance_unit = "m",
#'                    intern = TRUE)

dispersal_distance <- function(home_range = NULL,
                               home_range_unit = "ha",
                               dispersal_distance_unit = "m",
                               intern = TRUE){
  if(!is.numeric(home_range)){
    stop("Error. home_range must be numeric")
  }

  if(home_range_unit != "m2"){
    home_range <- unit_convert(home_range, home_range_unit, "m2")
  }

  dd <- data.frame("Type" = c("Median dispersal distance", "maximum dispersal distance"),
                   "Distance" = c(7*(sqrt(home_range)), 40*(sqrt(home_range))),
                   check.names = FALSE)

  if(dispersal_distance_unit != "m"){
    dd$Distance <- unit_convert(dd$Distance, "m", dispersal_distance_unit)
  }

  if(isTRUE(intern)){
    cat("Dispersal distance:\n"); print(dd)
  }

  return(dd)
}
