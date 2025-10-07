#' Estimation of connectivity indexes from CONEFOR
#'
#' Use the CONEFOR command line to estimate probabilistic and binary connectivity indexes
#' @param nodeFile \code{character}. Node file name. If the \code{prefix} parameter is used, then use the prefix name.
#' @param connectionFile \code{character}. Connection file name. If prefix is true use the prefix name.
#' @param folder  \code{character}. Path to folder workplace, default NULL.
#' @param coneforpath \code{character}. Path to \bold{Conefor 2.6 with command line interface}
#' (for more details and download \url{http://www.conefor.org/coneforsensinode.html}). Example, "C:/Users/coneforWin64.exe".
#' @param typeconnection \code{character}. Indicate if it is distance (dist), probability (prob), adjacency (adj).
#' @param typepairs \code{character}. \code{"all"} if all pairs of patches are in the distances file and \code{"notall"} if only they are the most important connections.
#' @param index \code{character}. Select Binary index (CCP, LCP, IIC, BC, BCIIC, NC, NL, H) or Probabilistic (F, AWF, PC, BCPC).
#' @param thdist \code{numeric}. Distance threshold to establish connections.
#' @param multdist \code{numeric}. Multiple distance threshold to establish connections (meters). The minimum (min) and maximum (max) distance values to be calculated need to be specified, as well as the value by which the distance value must be increased (step) starting from min until max is reached, e.g. \code{multdist = c(min=10000, max=100000, step=10000)}.
#' @param conprob \code{numeric}. If a probabilistic index is selected, indicate the probability of connection under the selected distance threshold, e.g., 0.5 that is 50 percentage of probability connection.
#' @param onlyoverall \code{logical}. Estimate only overall index, default \code{TRUE}.
#' @param LA \code{numeric}. Maximum landscape attribute.
#' @param nrestauration \code{logical}. If \code{TRUE} Indicates that there are nodes to add to the initial landscape -restauration effect-. The nodeFile needs a third column with values of 0 or 1 indicating for each node whether it already exists in the landscape (1) or is a candidate to be added to the landscape (0).The onlyoverall option have to be \code{FALSE}.
#' @param prefix \code{character}. Used in case of processing several sites at once.
#' @param write \code{character}. Full path direction of the folder.
#' @returns Four .txt files: "overall_indices", "results_all_EC(PC)", "results_all_overall_indices" y "node_importance".
#' @examples
#' \dontrun{
#' library(Makurhini)
#' setwd("~")
#' EstConefor(nodeFile = "nodes.txt",
#'           coneforpath = "C:/Users/coneforWin64.exe",
#'           connectionFile = "dist.txt",
#'           typeconnection = "dist",
#'           typepairs = "notall",
#'           index = "PC",
#'           thdist = 30000, multdist = NULL,
#'           conprob = 0.5,
#'           onlyoverall = FALSE, LA = NULL,
#'           nrestauration = FALSE,
#'           prefix= NULL, write = NULL)
#'}
#' @references
#' Saura, S. & Torné, J. 2012. Conefor 2.6 user manual (May 2012). Universidad Politécnica de Madrid.
#' Available at \url{www.conefor.org}.
#' @export
#' @importFrom utils read.table
EstConefor <- function(nodeFile,
                       coneforpath= NULL,
                       connectionFile,
                       folder = NULL,
                       typeconnection,
                       typepairs,
                       index,
                       thdist,
                       multdist = NULL,
                       conprob = NULL,
                       onlyoverall = TRUE,
                       LA = NULL,
                       nrestauration = FALSE,
                       prefix = NULL,
                       write = NULL){
  if(!file.exists(coneforpath)){
    stop("error, Conefor 2.6 with command line interface does not exist")
  }

  if(!file.exists(nodeFile)){
    stop("error, review nodeFile")
  }

  if(!file.exists(connectionFile)){
    stop("error, review connectionFile")
  }

  tt <- getwd()
  if(is.null(folder)){
    temp <- paste0(tempdir(), "\\TempCONEFOR", sample(1:1000, 1, replace = T))

    if(dir.exists(temp)){
      unlink(temp, recursive = TRUE)
    }

    dir.create(temp, recursive = TRUE)
    setwd(temp)
  } else {
    setwd(folder)
  }

  if (is.null(prefix)){
    p1 <- paste(conefor_exe, "-nodeFile", nodeFile, "-conFile", connectionFile, "-t", typeconnection, typepairs)
  } else {
    p1 <- paste(conefor_exe, "-nodeFile", nodeFile, "-conFile", connectionFile, "-t", typeconnection, typepairs, "-*")
  }

  binary <- c("CCP", "LCP", "IIC", "BC", "BCIIC", "NC", "NL", "H")
  probab <- c("F", "AWF", "PC", "BCPC")

  if (index %in% binary){
    if (index == "BCIIC"){
      p2 <- paste("-confAdj", thdist, paste0("-IIC -", index))
    }else {
      p2<-paste("-confAdj", thdist, paste0("-", index))
    }
  } else if (index %in% probab){
    if (index == "BCPC"){
      p2 <- paste("-confProb", thdist, conprob, paste0("-PC -", index))
    } else {
      p2 <- paste("-confProb", thdist, conprob, paste0("-", index))
    }
  } else {
    stop("error typeindex")
  }

  if (!is.null(multdist)){
    pM <- paste("-confMultipleValues", multdist[1], multdist[2], multdist[3], sep = " " )
  }

  if(isTRUE(onlyoverall) & isFALSE(nrestauration)){
    if (is.null(LA)){
      pO <- "onlyoverall"
    } else {
      pO <- paste("onlyoverall", "-landArea", LA)
    }

    if(is.null(prefix) & is.null(multdist)){
      result <- shell(paste(p1, p2, pO), intern = TRUE)
    } else if (is.null(prefix) & !is.null(multdist)){
      result <- shell(paste(p1, p2, pO, pM), intern = TRUE)
    } else if (!is.null(prefix) & is.null(multdist)){
      result <- shell(paste(p1, p2, pO, "-prefix", prefix), intern = TRUE)
    } else if (!is.null(prefix) & !is.null(multdist)){
      result <- shell(paste(p1, p2, pO, "-prefix", prefix, pM,  "-*", sep = " "), intern = TRUE)
    }
  } else if (isFALSE(onlyoverall) & isFALSE(nrestauration)){
    if (!is.null(LA)){
      pO <- paste("-landArea", LA)
    } else {
      pO <- ""
    }

    if (is.null(prefix) & is.null(multdist)){
      result <- shell(paste(p1, p2, pO), intern = TRUE)
    } else if (is.null(prefix) & !is.null(multdist)){
      result <- shell(paste(p1, p2, pO, pM), intern = TRUE)
    } else if (!is.null(prefix) & is.null(multdist)){
      result <- shell(paste(p1, p2, pO, "-prefix", prefix), intern = TRUE)
    } else if (!is.null(prefix) & !is.null(multdist)){
      result <- shell(paste(p1, p2, pO, "-prefix", prefix, pM,  "-*", sep = " "), intern = TRUE)
    }
  } else if (isFALSE(onlyoverall) & !isFALSE(nrestauration)){
    if (!is.null(LA)){
      pO <- paste("-landArea", LA)
    } else {
      pO<-""
    }
    if (is.null(prefix) & is.null(multdist)){
      result <- shell(paste(p1, p2, "-add", pO), intern = TRUE)
    } else if (is.null(prefix) & !is.null(multdist)){
      result <- shell(paste(p1, p2, "-add", pO, pM), intern = TRUE)
    } else if (!is.null(prefix) & is.null(multdist)){
      result <- shell(paste(p1, p2, "-add", pO, "-prefix", prefix), intern = TRUE)
    } else if (!is.null(prefix) & !is.null(multdist)){
      result <- shell(paste(p1, p2, "-add", pO, "-prefix", prefix, pM,  "-*", sep = " "), intern = TRUE)
    }
  }

  if (!isFALSE(unique(grepl("Error:", result, fixed = TRUE)))){
    setwd(tt)
    stop("Error. Please check the function parameters and input files")
  }

  #Export result
  df <- file.info(list.files(getwd(), full.names = TRUE))
  df <- df[order(df$atime, decreasing = TRUE),]

  df <- rownames(df)
  removefile_1 <- c(nodeFile, connectionFile, "coneforWin64.exe")

  for (i in 1:3){
    removefile <- which(basename(df) != removefile_1[[i]])
    df <- df[removefile]
  }

  result_2 <- tryCatch(lapply(1:length(df), function(x){
    x.1 <- df[x]
    x.1 <- read.table(x.1, header = if(x == 4){FALSE} else {TRUE})
    if(x == 4){
      names(x.1) <- c("Index", "Value")
    }
    return(x.1)
  }), error = function(err)err)

  if (inherits(result_2, "error")){
    removefile <- which(basename(df) != "overall_indices.txt")
    df <- df[removefile]
    result_2 <- tryCatch(lapply(1:length(df), function(x){
      x.1 <- df[x]
      x.1 <- read.table(x.1, header = if(x == 4){FALSE} else {TRUE})
      if(x == 4){
        names(x.1) <- c("Index", "Value")
      }
      return(x.1)
    }), error = function(err)err)
    nom <- basename(df)
    nom <- gsub(".txt","", nom)
    names(result_2) <- nom
  } else {
    nom <- basename(df)
    nom <- gsub(".txt","", nom)
    names(result_2) <- nom
  }

  if (!is.null(write)){
    setwd(tt)
    file.copy(df, write, overwrite = TRUE)
  }
  #Temove temporal files
  unlink(temp, recursive = TRUE)
  setwd(tt)
  return(result_2)
}


