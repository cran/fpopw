
###################################################################################################################################
#' @title select_Fpsn
#'
#' @description function to select the number of changepoints after Fpsn or Fpsn_w
#' using the penalty of Lebarbier 2005 given a estimator of the variance
#' @param res_fpsn output of Fpsn or Fpsn_w containg the costs in J.est and the segmented signal
#' @param method one of (1) "givenVariance" = using the penalty of Lebarbier 2005 given a estimator of the variance, (2) "biggest.S3IB" = biggest=TRUE in saut taken from S3IB, (3) "notbiggest.S3IB"  biggest=FALSE in saut taken from S3IB.
#' @param sigma variance used of the selection. If NULL use MAD on unweighted data.
#' @return return an integer: selected number of changes
#' @examples
#' x <- c(rnorm(100), rnorm(10^3)+2, rnorm(1000)+1)
#' res <- Fpsn_w(x=x, w=rep(1, length(x)), K=100)
#' select.res <- select_Fpsn(res, method="givenVariance")
#' smt <- getSMT(res, select.res)
#' @export
select_Fpsn <- function(res_fpsn, method="givenVariance", sigma=sdDiff(res_fpsn$signal)){
  K.seg <- 1:length(res_fpsn$J.est)
  n <- length(res_fpsn$signal)
  if(!is.null(res_fpsn$weights)) n <- sum(res_fpsn$weights > 0);
  
  pen.shape <- (2*K.seg*log(n/K.seg) + 5*(K.seg)) #*sigma^2

  if(!(method %in% c("givenVariance", "biggest.S3IB", "notbiggest.S3IB"))){
     warning("method should be on of givenVariance, biggest.S3IB, notbiggest.S3IB")
      return(NA)
  }
  ## strat sigma
  if(method == "givenVariance"){
    if(!is.numeric(sigma)){
      warning("a variance should be provided")
      return(NA)
    }
    K.sel <- which.min(res_fpsn$J.est + pen.shape*sigma^2) 
  }
  ## strat biggest
  if(method == "biggest.S3IB"){
    sel <- saut(-res_fpsn$J.est, pen.shape, K.seg, res_fpsn$n, TRUE)
    K.sel <- sel[1]
  }
  ## strat not-biggest
  if(method == "notbiggest.S3IB"){
    sel <- saut(-res_fpsn$J.est, pen.shape, K.seg, res_fpsn$n, TRUE)
    K.sel <- sel[1]
  }
  return(K.sel)
}

###################################################################################################################################
#' @title getTau_nomemory
#'
#' @description function to recover changes for a given selected K after fpsn_nomemory
#' @param res_fpsn output of the function res_fpsn_nomemory
#' @param K_selected K obtained using select_Fpsn
#' @return return a set of changes
#' @examples
#' x <- c(rnorm(100), rnorm(10^3)+2, rnorm(1000)+1)
#' res <- Fpsn_w_nomemory(x=x, w=rep(1, length(x)), K=100)
#' select.res <- select_Fpsn(res, method="givenVariance")
#' tau <- getTau_nomemory(res, select.res)
#' smt <- getSMT_(res$signal, res$weights, tau)
#' @export  
getTau_nomemory <- function(res_fpsn, K_selected){
  Kmax <- res_fpsn$Kmax
  K.seg <- 1:Kmax
  ## 
  ## recover segmentation with OP for a given penalty 
  if(K_selected == 1){
    tauHat <- length(res_fpsn$signal)
  }
  if(K_selected == Kmax){ ## largest
    ## TODO Warning 
    ## return with smallest penalty
    warning("selected K is equal to Kmax, select a smaller K")  
  }
  if(K_selected != 1 & K_selected != Kmax){
    ## return op with intermediate penalty
    betas <- (res_fpsn$J.est - res_fpsn$J.est[K_selected])/(K_selected - K.seg)
    K.smaller <- which.min(betas[1:(K_selected-1)])
    K.bigger  <- which.max(betas[(K_selected+1):Kmax]) + K_selected
    if(betas[K.smaller] <= betas[K.bigger]) { 
	## K_selected cannot be found with OP, should not happen 
	## if selection is done with a linear of concave penalty cannot happen
        ## in that case run with K.smaller
	## TODO warning
	warning("selected K cannot be recovered with OP (please select another K)")  
    } 
    beta.sel <- (res_fpsn$J.est[K.smaller] - res_fpsn$J.est[K.bigger])/(K.bigger - K.smaller)
    res.op <- fpopw::Fpop_w(res_fpsn$signal, res_fpsn$weights, beta.sel)
    tauHat <- res.op$t.est
  }
  return(tauHat)
}

###################################################################################################################################
#' @title saut
#'
#' @description model selection function taken from S3IB, 
#' @param Lv likelihood
#' @param pen penalty
#' @param Kseq number of changes
#' @param n number of datapoints
#' @param seuil threshold
#' @param biggest heuristic (biggest jump or slope)
#' @return a selected number of chagnes
#' @export
saut <- function(Lv, pen, Kseq, n, seuil=sqrt(n)/log(n), biggest=TRUE) {
   J=-Lv;Kmax=length(J); k=1;kv=c();dv=c();pv=c();dmax=1
   while (k<Kmax) {
	pk=(J[(k+1):Kmax]-J[k])/(pen[k]-pen[(k+1):Kmax])
	pm=max(pk); dm=which.max(pk); dv=c(dv,dm); kv=c(kv,k); pv=c(pv,pm)
	if (dm>dmax){  
	  dmax=dm; kmax=k; pmax=pm  
	}
	k=k+dm
  } 
  if (biggest){
	pv=c(pv,0); kv=c(kv,Kmax); dv=diff(kv); dmax=max(dv); rt=max(dv); rt=which(dv==rt)
	pmax=pv[rt[length(rt)]]
	alpha=2*pmax
	km=kv[alpha>=pv]; Kh =Kseq[km[1]] 
	return(c(Kh,alpha))
  } else {
	paux<-pv[which(kv<=seuil)]
	alpha<-2*min(paux)
	km=kv[alpha>=pv];	Kh =Kseq[km[1]] 
	return(c(Kh,alpha))	
  }
}


###################################################################################################################################
#' @title sdDiff
#'
#' @description Function to estimate the standard deviation
#' @param x signal
#' @param method used to estimate the variance : MAD or HALL
#' @return return a numeric value
#' @export
sdDiff <- function(x, method='MAD'){
  n = length(x)
  if(method == "MAD"){
    return(mad(diff(x)/sqrt(2)))	
  }
  if(method == "HALL"){
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(x)
    mat[2, -n] = mat[2, -1]
    mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]   
    return(sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3)))
  }
}

###################################################################################################################################
#' @title getSMT
#'
#' @description A function to get the smoothed profile from the output of Fpop, Fpop_w, Fpsn and Fpsn_w
#' @param res output of Fpop, Fpop_w, Fpsn or Fpsn_w
#' @param K the number of changes (only if Fpsn or Fpsn_w)
#' @return a vector of the smoothed profile
#' @export
getSMT <- function(res, K=NULL){
   ####################################################################################
   ### Warnings ...<
   if(!is.null(K) & res$method == "Fpop"){ ## if fpop and K -> don't use K
     warning("only one segmentation is available (fpop), parameter K will not be used")
   }
   if(is.null(K) & res$method == "Fpsn"){  ## if fpsn K should be given
     warning("Please provide a number of change to obtain the smt profile")
     return()
   }
   if(res$method == "FPSN_nomem"){  ## if fpsn_nomem DO NOT RUN FOR NOW
     warning("Please recover tau using tau <- getTau_nomemory(es_fpsn, K_selected)  and then apply getSMT_(x, weight, tau)")
     return()
   }
   ## >... Warnings
   ####################################################################################

   ### get the changepoints vector
   if(res$method == "Fpop"){ ## fpop
     pos <- res$t.est
   } else if(res$method == "Fpsn"){  ## fpsn
     pos <- res$t.est[K, 1:K]
   }

   ## recover the smt profil
   getSMT_(res$signal, res$weights, pos)
}

#' @title getSegSums_
#'
#' @description A function to get the segment sums of a vector given some changes including n
#' @param x data
#' @param tau changes (including n)
#' @return a vector of the sums
#' @export
getSegSums_ <- function(x, tau){
  csx <- cumsum(x)
  csx[tau] - c(0, csx[tau][-length(tau)])
}

#' @title getSMT_
#'
#' @description A function to get the smoothed profile from the data, weights and changepoints
#' @param x data
#' @param weights weights
#' @param tauHat changes (including n)
#' @return a vector of the smoothed profile
#' @export
getSMT_ <- function(x, weights=NULL, tauHat){
   if(is.null(weights)){
     weights <- rep(1, length(x))
   }

   lg.seg <- c(tauHat[1], diff(tauHat))
   sx <- getSegSums_(x*weights, tauHat)
   sw <- getSegSums_(weights, tauHat)
   smt    <- rep(sx/sw, lg.seg)
   return(smt)
}



