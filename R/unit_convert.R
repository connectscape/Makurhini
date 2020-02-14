#' Convert numeric types from one area unit to another
#'
#' @param data_unit numeric. Area value to be converted.
#' @param unit_1 character. Select an area unit: "m2" = square meter, "Dam2" = square decametre,
#' "ha" = hectares, "km2" = square kilometer, "inch2" = square inch, "foot2" = square foot,
#' "yard2" = square yard, "mile2" = square mile. Length unit: "m" = meter, "km" = kilometer, "inch", "foot", "yard", "mile"
#' @param unit_2 character. Select an area unit: "m2" = square meter, "Dam2" = square decametre,
#' "ha" = hectares, "km2" = square kilometer, "inch2" = square inch, "foot2" = square foot,
#' "yard2" = square yard, "mile2" = square mile. Length unit: "m" = meter, "km" = kilometer, "inch", "foot", "yard", "mile"
#' @export

unit_convert <- function(data_unit,
                         unit_1 = c("m2", "Dam2", "ha", "km2", "inch2", "foot2", "yard2", "mile2", "m", "km", "inch", "foot", "yard", "mile"),
                         unit_2 = c("m2", "Dam2", "ha", "km2", "inch2", "foot2", "yard2", "mile2", "m", "km", "inch", "foot", "yard", "mile")){

  if(unit_1 %in% c("m2", "Dam2", "ha", "km2", "inch2", "foot2", "yard2", "mile2") &
     unit_2 %in% c("m2", "Dam2", "ha", "km2", "inch2", "foot2", "yard2", "mile2")){
    if(unit_1== "m2"){
      uconv <- data.frame(m2 = 1, km2 =  0.000001, ha = 0.0001, Dam2 = 0.01, inch2 = 1550, foot2 =10.7639, yard2 = 1.19599, mile2 = 3.861e-7)
    } else if (unit_1 == "Dam2"){
      uconv <- data.frame(km2 = 0.0001, ha = 0.01, Dam2 = 1, m2 = 100, inch2 = 155000, foot2 =1076.39, yard2 = 119.599, mile2 = 3.861e-5)
    } else if(unit_1 == "ha"){
      uconv <- data.frame(km2 = 0.01, ha = 1, Dam2 = 100, m2 = 10000, inch2 = 1.55e+7, foot2 = 107639, yard2 = 11959.9, mile2 = 0.00386102)
    } else if(unit_1 == "km2"){
      uconv <- data.frame(km2 = 1, ha = 100, Dam2 = 10000, m2 = 1000000, inch2 = 1.55e+9, foot2 = 1.076e+7,
                          yard2 = 1.196e+6, mile2 = 0.386102)
    } else if(unit_1 == "inch2"){
      uconv <- data.frame(km2 = 6.4516e-10, ha = 6.4516e-8, Dam2 = 6.4516e-6,  m2 = 0.00064516, inch2 = 1,
                          foot2 = 0.00694444, yard2 = 0.000771605, mile2 = 2.491e-10)
    } else if(unit_1 == "foot2"){
      uconv <- data.frame(km2 = 9.2903e-8, ha = 9.2903e-6, Dam2 = 0.00092903, m2 = 0.092903,
                          inch2 = 144, foot2 = 1, yard2 = 0.111111, mile2 = 3.587e-8)
    } else if(unit_1 == "yard2"){
      uconv <- data.frame(km2 = 8.3613e-7, ha = 8.3613e-5, Dam2 = 0.00836127, m2 = 0.836127,
                          inch2 = 1296, foot2 = 9, yard2 = 1, mile2 = 3.2283e-7)
    } else {
      uconv <- data.frame(km2 = 2.58999, ha = 258.999, Dam2 = 25899.9,  m2 = 2.59e+6,
                          inch2 = 4.014e+9, foot2 = 2.788e+7, yard2 = 3.098e+6, mile2 = 1)
    }
  } else if (unit_1 %in% c("m", "km", "inch", "foot", "yard", "mile") &
             unit_2 %in% c("m", "km", "inch", "foot", "yard", "mile")){
    if(unit_1== "m"){
      uconv <- data.frame(m = 1, km = 0.001, inch = 39.3701, foot = 3.28084, yard = 1.09361, mile = 0.000621371)
    } else if (unit_1 == "km"){
      uconv <- data.frame(m = 1000, km = 1, inch = 39370.1, foot = 3280.84, yard = 1093.61, mile = 0.621371)
    } else if(unit_1 == "inch"){
      uconv <- data.frame(m = 0.0254, km = 2.54e-5, inch = 1, foot = 0.0833333, yard = 0.0277778, mile = 1.578282197e-5)
    } else if(unit_1 == "foot"){
      uconv <- data.frame(m =  0.3048, km = 0.0003048, inch = 12, foot = 1, yard = 0.333333, mile = 0.000189394)
    } else if(unit_1 == "yard"){
      uconv <- data.frame(m = 0.9144, km = 0.0009144, inch = 36, foot = 3, yard = 1, mile =0.000568182)
    } else {
      uconv <- data.frame(m = 1609.34, km = 1.60934, inch = 63360, foot = 5280, yard = 1760, mile = 1)
      }
  } else {
    stop("Select a correct unit area")
  }

  x1 <- data_unit * uconv[[unit_2]]

  return(x1)
}

