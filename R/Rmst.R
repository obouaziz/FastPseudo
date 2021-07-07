#' Restricted Mean Survival Time
#'
#' Compute the Restricted Mean Survival Time (RMST) for right-censored data
#'
#' @param Time a continuous time variable.
#' @param status an indicator for censoring for the corresponding Time variable.
#' @param tau the time endpoint for computing the RMST

#' @export
Rmst <- function(Time,status,tau) {UseMethod("Rmst")}
#' @export
Rmst.default<-function(Time,status,tau){
  if (length(Time)!=length(status)) stop("'Time' and 'status' must have the same length")
  if (is.null(tau)) stop("a finite value of tau must be specified")
  if (is.na(tau)|(tau==Inf)) stop("a finite value of tau must be specified")
  if (tau<0) stop("'tau' must be positive")
  Tsort<-sort(Time,index.return=TRUE)
  Tobs_ord<-Tsort$x
  status_ord<-status[Tsort$ix]
  result<-list(rmst=Rmst_ord(Tobs_ord,status_ord,tau),tau=tau)
  class(result)<-"Rmst"
  return(result)
}
#' @export
print.Rmst <- function(x, ...)
{
  cat(paste(" restricted mean with upper limit = ",x$tau,"\n",sep=""))
  print(round(x$rmst[1],4))#x$rmst[1]
}

### Helper ####
Rmst_ord=function(Tobs_ord,status_ord,tau){
  resKM=survival::survfit(survival::Surv(Tobs_ord,status_ord)~1)
  #for large sample sizes, sometimes some Tobs_ord are removed in survfit...
  #then length(resKM$surv) is <n and the codes won't work
  #this is a fix
  if (length(resKM$surv)!=length(Tobs_ord)){
    KMfun=stats::stepfun(resKM$time,c(1,resKM$surv))
    resKMsurv=KMfun(Tobs_ord)
  } else {resKMsurv=resKM$surv}
  if (max(Tobs_ord[status_ord])<tau){ #case where tau is greater than last observed event of interest
    K=sum(status_ord)
    integ=c(resKMsurv[status_ord][1:(K-1)])*diff(Tobs_ord[status_ord])
    return(Tobs_ord[status_ord][1]+sum(integ)+resKMsurv[status_ord][K]*(tau-Tobs_ord[status_ord][K]))
  } else {
    #case where tau is less than last observed event of interest
    if (Tobs_ord[status_ord][1]>=tau) { #case where tau is less than first observed event of interest
      return(tau)
    } else{#other cases
      index=which(Tobs_ord[status_ord]>=tau)[1]-1 #separate integral in three parts
      integ=c(resKMsurv[status_ord])[1:(index-1)]*diff(Tobs_ord[status_ord])[1:(index-1)] #integral between T[1,delta] and last observed event before tau
      return(Tobs_ord[status_ord][1]+sum(integ)+(tau-Tobs_ord[status_ord][index])*c(resKMsurv[status_ord])[index])
      #the first part is for the integral between 0 and T[1,delta]
      #the last part is for the integral between last observed event before tau and tau
    }
  }
}
