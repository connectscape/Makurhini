#' Multiple ProtConn statistic
#'
#' @param ProtConn ProtConn table (see, ProtConnCLaMult())
#' @param ci character. A character vector representing the type of confidence intervals that will be estimated. The value should be any subset of the values c("norm","basic", "stud", "perc", "bca") or "all" which will compute all five types of intervals (see, boot::boot.ci())
#' @param nr numeric. The number of bootstrap replicates
#' @references Carpenter, J., & Bithell, J. (2000). Bootstrap con " dence intervals : when , which , what ? A practical guide for medical statisticians. Statist. Med., 19, 1141–1164.
#' @export
#' @importFrom boot boot boot.ci
#' @importFrom purrr map
#' @importFrom stats sd
ProtConnStat <- function(ProtConn, ci= "all", nr=500){
  if(!is.null(ci)){
    if(ci == "all"){
      ci= c("normal", "basic", "percent", "bca")
    } }

  ProtConn <- ProtConn[c(which(names(ProtConn) == "EC(PC)"):ncol(ProtConn))]
  ProtConn[,1] <- as.character(ProtConn[,1]) %>% as.numeric()
  ProtConn[,2] <- as.character(ProtConn[,2]) %>% as.numeric()

  ProtConn[is.na(ProtConn)] = 0
  ProtConn2 <- cbind(colMeans(ProtConn))
  ProtConn2 <- as.data.frame(cbind(ProtConn2, cbind(apply(ProtConn,2, function(x) sd(x))),
                                   cbind(apply(ProtConn,2, function(x) sd(x)/sqrt(length(x))))))
  colnames(ProtConn2)<-c("Values", "SD", "SEM")
  mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)

  #IC
  if(!is.null(ci)){
    ICBoot <- map(as.list(ProtConn) , function(n){
      if(length(unique(n)) > 1){
        set.seed(626)
        bootcorr <- boot(n, mean.fun, R = nr)
        err2 <- tryCatch(boot.ci(boot.out = bootcorr, type = c("norm", "basic", "perc", "bca")),
                         error = function(err)err)
        if (inherits(err2, "error")){
          err2 <- rep(n, length(ci) * 2)
        } else {
          err2 <-err2[ci]
          err2 <-lapply(err2, function(x) as.data.frame(x)[(length(x)-1):length(x)])
          for (j in 1:length(err2)){
            colnames(err2[[j]]) <- c("lower", "upper")
          }
          err2 <-do.call(cbind, err2)
        }
      } else {
        err2 <- rep(unique(n), 8)
        names(err2) <- c("normal.lower", "normal.upper", "basic.lower", "basic.upper",
                         "percent.lower", "percent.upper", "bca.lower", "bca.upper")
        err2 <- map(as.list(ci), function(x){
          x2 <- err2[which(startsWith(names(err2), x))]
          x2 <- as.data.frame(x2)
          return(x2)})
        err2 <- do.call(rbind, err2) %>% t()
        rownames(err2) <- 1
      }
      return(err2)})
    ICBoot <- do.call(rbind, ICBoot)
    ICBoot <- apply(ICBoot, 2, function(x){
      x[which(x >100)]<- 100
      x[which(x <0)]<- 0
      return(x)})
    ProtConn2 <- cbind(ProtConn2, ICBoot)
  }
  return(ProtConn2)
}

