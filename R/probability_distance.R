#'Plot probability  of dispersal
#'
#' Negative exponential dispersal kernel to calculate the probability of dispersal between two nodes. Used in probabilistic indexes (e.g., PC, BCPC, ProtConn etc.)
#' @param probability numeric. Probability of dispersal at median_distance. If NULL
#' @param median_distance numeric. Up to five maximum dispersal distance (km)
#' @param eval_distance numeric. Calculate the probability of dispersal at a specific median distances (km). Available when only one median_distance is used.
#' @param min.prob numeric. Value between 0-1, maximum x axe value
#' @examples
#' \dontrun{
#' probability_distance(probability= 0.5, median_distance = c(1, 10, 30, 100), min.prob = 0.01)
#' probability_distance(probability= 0.5, median_distance = 100, eval_distance = 100)
#' probability_distance(probability= 0, median_distance = 100, eval_distance = 1)
#' }
#' @export
#' @importFrom purrr map map_dbl
#' @importFrom graphics par plot axis box lines legend

probability_distance <- function(probability, median_distance, eval_distance = NULL,
                                 min.prob = 0.01){

  if(probability == 0){
    probability = 1e-100
  }
  x <-1
  repeat{
    x = x+1
    m <- exp(x * log(probability)/max(median_distance))
    if (m < min.prob & x > max(median_distance)){
      break
      return(x)
    }
  }

  m = x
  Rcolors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

  if(!is.null(eval_distance)){
    if(m < eval_distance){
      m <- eval_distance + eval_distance/100
    }

    if(length(median_distance) == 1){
      data1 <- exp(1:m * log(probability)/median_distance)

      result_1 <- exp(eval_distance * log(probability)/median_distance)

      t <- c(rep(result_1, eval_distance - 1), result_1, 0)

      par(xaxs="i")
      plot(data1, type = "l", xlab= paste0("Internode distance", " (dij)"), ylab = "Probability of dispersal (pij)",
           lwd = 2, xlim=c(0,m), axes = F, main = paste0("pij,", eval_distance, " = ", round(result_1, 10)))
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()
      lines(t, type = "l", col = Rcolors[1], lty = 2, lwd = 2)
      legend("topright", c(paste0("Distance ", median_distance), "Evaluated distance"), lwd=c(2,2), lty = c(1,2), col=c("black",Rcolors[1]), y.intersp=1.5)
      return(result_1)
    } else {
      data_list <- map(as.list(median_distance), function(x){
        data1 <- exp(1:m * log(probability)/x)
        result_1 <- exp(eval_distance * log(probability)/x)
        t <- c(rep(result_1, eval_distance - 1), result_1, 0)
        return(list(data1, result_1, t))})

      maxv <- map_dbl(data_list, max) %>% max()
      minv <- map_dbl(data_list, min) %>% min()

      par(xaxs="i")
      plot(data_list[[1]][[1]], type = "l", xlab= paste0("Internode distance ", "(dij)"),
           ylab = "Probability of dispersal (pij)",
           ylim = c(minv, maxv),
           lwd = 2, xlim=c(0,m), axes=F)
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()
      lines(data_list[[1]][[3]], type = "l", col = Rcolors[1], lty = 2, lwd = 2)

      for(i in 2:length(median_distance)){
        lines(data_list[[i]][[1]], type = "l", col = Rcolors[i], lwd = 2)
        lines(data_list[[i]][[3]], type = "l", col = Rcolors[1], lty = 2, lwd = 2)
      }
      legend("topright", c(paste0("Distance ", median_distance), "Evaluated distance"), lwd=2, lty = c(rep(1,length(median_distance)), 2),
             col= c("black", Rcolors[2:length(median_distance)], Rcolors[1]), y.intersp=1.5, cex=0.8)

      result_1 <- map(data_list, function(x){x[[2]]})
      names(result_1) <- paste0("d_", median_distance)
      return(result_1)
    }
  } else {
    if(length(median_distance) == 1){
      data <- exp(1:m * log(probability)/median_distance)
      par(xaxs="i")
      plot(data, type = "l", xlab=  "Distance (dij)", ylab = "Probability of dispersal (pij)",  axes=F)
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()
      lines(data_list[[1]][[3]], type = "l", col = Rcolors[1], lty = 2, lwd = 2)
      legend("topright", paste0("Distance ", median_distance), lwd=c(2,2), lty = c(1,2), col=c("black",Rcolors[1]), y.intersp=1.5)
    } else {
      data_list <- map(as.list(median_distance), function(x){
        data1 <- exp(1:m * log(probability)/x)
        return(data1)})

      maxv <- map_dbl(data_list, max) %>% max()
      minv <- map_dbl(data_list, min) %>% min()

      par(xaxs="i")
      plot(data_list[[1]], type = "l", xlab= "Distance (dij)",
           ylab = "Probability of dispersal (pij)",
           ylim = c(minv, maxv),
           lwd = 2, xlim=c(0,m), axes=F)
      axis(side = 1, at= round(seq(0,m, m/10)))
      axis(side = 2, at= seq(0, 1, 1/10))
      box()

      for(i in 2:length(median_distance)){
        lines(data_list[[i]], type = "l", col = Rcolors[i-1], lwd = 2)
      }
      legend("topright", paste0("Distance ", median_distance), lwd=2,
             col= c("black", Rcolors[1:length(median_distance)]), y.intersp=1.5, cex = 0.8)
    }
  }

}
