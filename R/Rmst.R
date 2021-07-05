#' Restricted Mean Survival Time
#'
#' Compute the Restricted Mean Survival Time (RMST) for right-censored data
#'
#' @param Time a continuous time variable. Must be sorted
#' @param status an indicator for censoring for the corresponding Time variable.
#' @param tau the time endpoint for computing the RMST

Rmst=function(Time,status,tau){
  Tobs_ord<-Time
  status_ord<-status
  resKM=survival::survfit(survival::Surv(Tobs_ord,status_ord)~1)
  if (max(Tobs_ord[status_ord])<tau){ #case where tau is greater than last observed event of interest
    K=sum(status_ord)
    integ=c(resKM$surv[status_ord][1:(K-1)])*diff(Tobs_ord[status_ord])
    return(Tobs_ord[status_ord][1]+sum(integ)+resKM$surv[status_ord][K]*(tau-Tobs_ord[status_ord][K]))
  } else {
    #case where tau is less than last observed event of interest
    if (Tobs_ord[status_ord][1]>=tau) { #case where tau is less than first observed event of interest
      return(tau)
    } else{#other cases
      index=which(Tobs_ord[status_ord]>=tau)[1]-1 #separate integral in three parts
      integ=c(resKM$surv[status_ord])[1:(index-1)]*diff(Tobs_ord[status_ord])[1:(index-1)] #integral between T[1,delta] and last observed event before tau
      return(Tobs_ord[status_ord][1]+sum(integ)+(tau-Tobs_ord[status_ord][index])*c(resKM$surv[status_ord])[index])
      #the first part is for the integral between 0 and T[1,delta]
      #the last part is for the integral between last observed event before tau and tau
    }
  }
}
