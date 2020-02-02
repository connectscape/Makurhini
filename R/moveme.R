#' Move columns
#'
#' @param invec A matrix of values with at least two dimensions e.g. data frame, attributes of a shapefile
#' @param movecommand character. Instructions to move the column, four options:first, last, before and after. e.g.
#' @examples
#' \dontrun{ moveme(names(df), "g first")
#'  moveme(names(df), "g first; a last; e before c")
#'  df[moveme(names(df), "g first")]
#'  }
#' @references
#' The solution was proposed by the user "A5C1D2H2I1M1N2O1R2T1" \url{https://stackoverflow.com/users/1270695/a5c1d2h2i1m1n2o1r2t1}.
#' @export
moveme <- function (invec, movecommand){
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]],",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x){
    Where <- x[which(x %in% c("before", "after", "first", "last")):length(x)]
    ToMove <- base::setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- base::setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      } else if (A == "after") {
        after <- match(ba, temp)
      }
    } else if (A == "first") {
      after <- 0
    } else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}
