#' Function to contrast differences of dA and dECA
#' @param dECA dECA
#' @param dA dA

dECAfun<-function(dECA, dA) {
  if (dECA < dA & dA < 0){
    TC <- "dECA < dA < 0"
    } else if (dA < dECA & dECA < 0){
      TC <- "dA < dECA < 0"
      } else if (dECA == dA & dA < 0) {
        TC <- "dECA = dA < 0"
        } else if (dECA == dA & dA == 0) {
          TC <- "dECA = dA = 0"
          }else {
          TC <- "dECA or dA gain"
        }
  return(TC)
}

#' Function to contrast differences of dA and dECA
#' @param dECA dECA
#' @param dA dA

dECAfun2<-function(dECA, dA) {
  if (dECA < dA & dA < 0){
    TC <- "+ Connectivity loss"
} else if (dA < dECA & dECA < 0){
    TC <- "+ Habitat loss"
  } else if (dECA == dA & dA < 0) {
    TC <- "Equal loss"
  } else if (dECA == dA & dA == 0) {
    TC <- "Habitat or connectivity maintained"
  } else {
    TC <- "Habitat or connectivity gain"
    }
  return(TC)
}
