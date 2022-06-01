###################################################################################################################################
#' @title retour_op
#'
#' @description Function used internally by Fpop and Fpop_w to do the backtracking and recover the best set of changes from 1 to i 
#' @param path vector of length n containing the best last changes for any j in \eqn{[1, n]}. This vector is computed in the Fpop and Fpop_w using the colibri_op_c or colibri_op_weight_c function.
#' @param i the last position to consider to start the backtracking.
#' @return set of optimal changes up to i.
#' @export
retour_op <- function(path, i){
   chaine <- integer(1)
   chaine[1] <- length(path)
   j <- 2
   while(chaine[j-1] > 0){
	chaine[j] <- path[chaine[j-1]]
	j=j+1	
	}
   return(rev(chaine)[-1]) 
}

###################################################################################################################################
#' @title retour_sn
#'
#' @description Function used internally by Fpsn and Fpsn_w to do the backtracking and recover the best set of segmentations in 1 to K changes from 1 to n.
#' @param path matrix of size (K x n) containing the last optimal changes up to j in k segments with i in \eqn{[1, n]} and k in \eqn{[1, K]}. This matrix is computed in the Fpsn or Fpsn_w function using the colibri_sn_c or colibri_sn_weight_c functions.
#' @return a matrix of size (K x K) containing the best segmentations in 1 to K segments.
#' @export
retour_sn <- function(path){
  n <- ncol(path)
  res3 <- matrix(NA, nrow=nrow(path), ncol=nrow(path))
  res3[1, 1] <- 0
  for(i in 2: nrow(path)){
    res3[i, i-1] <- path[i, n]
    for(k in 1:(i-1)){
      res3[i, i-1-k] <- path[i-k, res3[i, i-k]]
    }
    
  }
  diag(res3) <- ncol(path)
  return(res3)
}

###################################################################################################################################
#' @title Fpop
#'
#' @description Function to run the Fpop algorithm (Maidstone et al. 2016). It uses functional pruning and optimal partionning. It optimizes the L2-loss for a penalty lambda per change.
#' @param x a numerical vector to segment
#' @param lambda the penalty per changepoint (see Maidstone et al. 2016) 
#' @param mini minimum mean segment value to consider in the optimisation.
#' @param maxi maximum mean segment value to consider in the optimisation.
#' @return return a list with a vector t.est containing the position of the change-points, the number of changes K and, the cost J.est.
#' @examples 
#' x <- c(rnorm(100), rnorm(10^3)+2, rnorm(1000)+1)
#' est.sd <- sdDiff(x) ## rough estimate of std-deviation
#' res <- Fpop(x=x,lambda=2*est.sd^2*log(length(x)))
#' smt <- getSMT(res)
#' @export
Fpop <- function(x, lambda, mini=min(x), maxi=max(x)){
  n <- length(x)
  A <- .C("colibri_op_R_c", signal=as.double(x), n=as.integer(n), 
		lambda=as.double(lambda),   min=as.double(mini), 
		max=as.double(maxi), path=integer(n), cost=double(n)
	, PACKAGE="fpopw")
    A$t.est <- retour_op(A$path, n)
    A$K <- length(A$t.est)
    A$J.est <- A$cost[n] - (A$K+1)*lambda + sum((x^2))
    A$weights <- NULL ## used to compute the smt profile
    A$method <- "Fpop"
    return(A);	
} 

###################################################################################################################################
#' @title Fpop_w
#'
#' @description Function to run the Fpop algorithm (Maidstone et al. 2016) with weights. It uses functional pruning and optimal partionning. It optimizes the weighted L2-loss (\eqn{w_i (x_i - \mu)2}) for a penalty lambda per change.
#' @param x a numerical vector to segment.
#' @param w a numerical vector of weights (values should be larger than 0).
#' @param lambda the penalty per changepoint (see Maidstone et al. 2016).
#' @param mini minimum mean segment value to consider in the optimisation.
#' @param maxi maximum mean segment value to consider in the optimisation.
#' @return return a list with a vector t.est containing the position of the change-points, the number of changes K and, the cost J.est.
#' @examples 
#' x <- c(rnorm(100), rnorm(10^3)+2, rnorm(1000)+1)
#' est.sd <- sdDiff(x) ## rough estimate of std-deviation
#' res <- Fpop_w(x=x, w=rep(1, length(x)), lambda=2*est.sd^2*log(length(x)))
#' smt <- getSMT(res)
#' @export
Fpop_w <- function(x, w, lambda, mini=min(x), maxi=max(x)){
  n <- length(x)
  if(min(w) < 0){
    warning("All weights should be larger than 0")
    return()
  }
  A <- .C("colibri_op_weights_R_c", signal=as.double(x), weights=as.double(w), n=as.integer(n), 
		lambda=as.double(lambda),   min=as.double(mini), 
		max=as.double(maxi), path=integer(n), cost=double(n)
	, PACKAGE="fpopw")
  A$t.est <- retour_op(A$path, n)
  A$K <- length(A$t.est)
  A$J.est <- A$cost[n] - (A$K+1)*lambda + sum(w*(x^2))
  A$method <- "Fpop"
  return(A);	
} 



###################################################################################################################################
#' @title Fpsn
#'
#' @description Function to run the pDPA algorithm (Rigaill 2010 and 2015). It uses functional pruning and segment neighborhood. It optimizes the L2-loss for 1 to Kmax changes.
#' @param x a numerical vector to segment
#' @param Kmax max number of segments (segmentations in 1 to Kmax segments are recovered).
#' @param mini minimum mean segment value to consider in the optimisation
#' @param maxi maximum mean segment value to consider in the optimisation
#' @return return a list with a matrix t.est containing the change-points of the segmentations in 1 to Kmax changes and, the cost J.est in 1 to Kmax changes.
#' @examples 
#' x <- c(rnorm(100), rnorm(10^3)+2, rnorm(1000)+1)
#' res <- Fpsn(x=x, K=100)
#' select.res <- select_Fpsn(res, method="givenVariance")
#' smt <- getSMT(res, select.res)
#' @export
Fpsn <- function(x, Kmax, mini=min(x), maxi=max(x)){
  n <- length(x)
  A <- .C("colibri_sn_R_c", signal=as.double(x), n=as.integer(n), 
		Kmax=as.integer(Kmax),   min=as.double(mini), 
		max=as.double(maxi), path=integer(Kmax*n), J.est=double(Kmax), allCost=double(n*Kmax)
	, PACKAGE="fpopw")
    A$path <- matrix(A$path, nrow=Kmax, byrow=TRUE)
    A$allCost <- matrix(A$allCost, nrow=Kmax, byrow=TRUE)
    A$t.est <- retour_sn(A$path)
    A$weights <- NULL ## used to compute the smt profile
    A$method <- "Fpsn"
    return(A);	
} 



###################################################################################################################################
#' @title Fpsn_w
#'
#' @description Function to run the weighted pDPA algorithm (Rigaill 2010 and 2015). It uses functional pruning and segment neighborhood. It optimizes the weighted L2-loss (\eqn{w_i (x_i - \mu)2}) for 1 to Kmax changes.
#' @param x a numerical vector to segment
#' @param w a numerical vector of weights (values should be larger than 0).
#' @param Kmax max number of segments (segmentations in 1 to Kmax segments are recovered).
#' @param mini minimum mean segment value to consider in the optimisation
#' @param maxi maximum mean segment value to consider in the optimisation
#' @return return a list with a matrix t.est containing the change-points of the segmentations in 1 to Kmax changes and, the costs J.est in 1 to Kmax changes.
#' @examples 
#' x <- c(rnorm(100), rnorm(10^3)+2, rnorm(1000)+1)
#' res <- Fpsn_w(x=x, w=rep(1, length(x)), K=100)
#' select.res <- select_Fpsn(res, method="givenVariance")
#' smt <- getSMT(res, select.res)
#' @export
Fpsn_w <- function(x, w, Kmax, mini=min(x), maxi=max(x)){
  n <- length(x)
  if(min(w) < 0){
    warning("All weights should be larger than 0")
    return()
  }
  A <- .C("colibri_sn_weights_R_c", signal=as.double(x), weights=as.double(w), n=as.integer(n), 
		Kmax=as.integer(Kmax),   min=as.double(mini), 
		max=as.double(maxi), path=integer(Kmax*n), J.est=double(Kmax), allCost=double(n*Kmax)
	, PACKAGE="fpopw")
  A$path <- matrix(A$path, nrow=Kmax, byrow=TRUE)
  A$t.est <- retour_sn(A$path)
  A$allCost <- matrix(A$allCost, nrow=Kmax, byrow=TRUE)
  A$method <- "Fpsn"
  return(A);	
} 




###################################################################################################################################
#' @title Fpsn_w_nomemory
#'
#' @description Function to run the weighted pDPA algorithm (Rigaill 2010 and 2015) without storing the set of last changes. It only return the cost in 1 to Kmax changes. It uses functional pruning and segment neighborhood. It optimizes the weighted L2-loss (\eqn{w_i (x_i - \mu)2}) for 1 to Kmax changes.
#' @param x a numerical vector to segment
#' @param w a numerical vector of weights (values should be larger than 0).
#' @param Kmax max number of segments (segmentations in 1 to Kmax segments are recovered).
#' @param mini minimum mean segment value to consider in the optimisation
#' @param maxi maximum mean segment value to consider in the optimisation
#' @return return a list with the costs J.est in 1 to Kmax changes.
#' @examples 
#' res <- Fpsn_w_nomemory(x=rnorm(10^4), w=rep(1, 10^4), K=100)
#' @export
Fpsn_w_nomemory <- function(x, w, Kmax, mini=min(x), maxi=max(x)){
  n <- length(x)
  A <- .C("colibri_sn_weights_nomemory_R_c", signal=as.double(x), weights=as.double(w), n=as.integer(n), 
		Kmax=as.integer(Kmax),   min=as.double(mini), 
		max=as.double(maxi), J.est=double(Kmax)
	, PACKAGE="fpopw")
  A$method <- "Fpsn_nomem"
  return(A);	
} 


