#' Close cores
#'
#' @param w numeric. workers or number of cores
#' @references Revanth Nemani. 2020. https://stackoverflow.com/questions/59652303/kill-futures-from-future-apply-on-linux#comment105470849_59652303
#' @importFrom tools pskill
#' @importFrom ps ps
#' @export
close_multiprocess <- function(w){
  w <- as.numeric(w)
  a <- ps::ps()
  a <- a[which(a$name =="Rscript.exe"),1]

  for(i in a){
    tools::pskill(i)
  }
}
