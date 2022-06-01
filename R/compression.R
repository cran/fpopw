
#' @title compress.data
#'
#' @description compress data and return a weighted profile 
#' @param x a numerical vector
#' @return a list with the compressed profile x and associated repeat vector vrep
#' @export
compress.data <- function(x){
  n <- length(x)
  tau.CP <- c(which(diff(x) != 0), n)
  x.CP   <- x[tau.CP]
  vec.rep   <- diff(c(0, tau.CP))
  return(list(x=x.CP, vec.rep=vec.rep))
}

#' @title uncompress.vec
#'
#' @description return a vector to uncompress a profile, segmentation or smt
#' @param vec.rep integer vector with the number of time each point should be repeated 
#' @return return a vector to uncompress a profile, segmentation or smt
#' @export
uncompress.vec <- function(vec.rep){
  rep(1:length(vec.rep), vec.rep)
}

#' @title decompress.smt
#'
#' @description vector to decompress a compressed smoothed profile (a call to rep)
#' @param smt.CP smoothed and compressed profile
#' @param vec.rep weights to use for decompression
#' @return a vector to replicate duplicated datapoints
#' @export
uncompress.smt <- function(smt.CP, vec.rep){
  rep(smt.CP, vec.rep)
}

#' @title get.change
#'
#' @description Function returning changes in a smoothed profile
#' @param smt smoothed profile
#' @return a vector of changes including n
#' @export
get.change <- function(smt){
  c(which(diff(smt) !=0), length(smt))
}

