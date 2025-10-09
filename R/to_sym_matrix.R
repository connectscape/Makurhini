#' Pairwise to symmetric matrix
#' @param df data.frame
#' @param from character
#' @param to character
#' @param value character
#' @param diag_value numeric
#' @param fill_na numeric
#' @keywords internal

to_sym_matrix <- function(df, from = "from", to = "to", value = "value", diag_value = 0, fill_na = 0) {
  f <- df[[from]]; t <- df[[to]]; v <- df[[value]]
  nodes <- sort(unique(c(f, t)))
  m <- matrix(NA_real_, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  diag(m) <- diag_value
  for (k in seq_len(nrow(df))) {
    i <- match(f[k], nodes); j <- match(t[k], nodes); val <- v[k]
    m[i, j] <- val; m[j, i] <- val
  }
  m[is.na(m)] <- fill_na
  return(m)
}
