#' Restricted Mean Survival Time for Interval-Censored data
#'
#' Compute the Restricted Mean Survival Time (RMST) for interval-censored data
#' from the piecewise constant hazard (PCH) model
#'
#' @param cuts the cuts used in the PCH model.
#' @param alpha the value of the hazard function between each cut. Should be of
#' length equal to \code{length(cuts)}+1.
#' @param tau the time endpoint for computing the RMST.
#' @return returns the RMST along with the value \code{tau} that was used.
#' @seealso \code{\link{pseudoIC}}, \code{\link{mleIC}}, \code{\link{pchcumhaz}}, \code{\link{Rmst}}.

#' @export
RmstIC <- function(cuts,alpha,tau) {UseMethod("RmstIC")}
#' @export
RmstIC.default<-function(cuts,alpha,tau){
  if (length(alpha)!=(length(cuts)+1)) stop("'alpha' must be one length longer than 'cuts'")
  if (is.null(tau)) stop("a finite value of tau must be specified")
  if (is.na(tau)|(tau==Inf)) stop("a finite value of tau must be specified")
  if (tau<0) stop("'tau' must be positive")
  a=log(alpha)
  if (is.null(cuts))
  {
    result=exp(log1p(-exp(-alpha*tau))-a)
  } else {
    index=c(tau>cuts) #watch out that tau>0!c(TRUE,tau>cuts)
    nb_index=sum(index)
    k=length(a)
    cuts0=c(0,cuts)
    if (nb_index==0){
      cuts_tau=c(0,tau)
    } else {
      cuts_tau=c(0,cuts[index],tau)#rep(tau,k-nb_index-1)#Note that length(cuts_tau)=nb_index+2
    }
    part1=-c(0,cumsum(alpha*c(diff(cuts0),0))[-k])#+cuts0*exp(a)
    part1=part1[1:(nb_index+1)]
    part2=log1p(-exp(-alpha[1:(nb_index+1)]*diff(cuts_tau)))-a[1:(nb_index+1)]
    result=(exp(part1+part2))
    if (sum(a[1:(nb_index+1)]==-Inf)>0) {
      aInf<-(a[1:(nb_index+1)]==-Inf)
      result[aInf]<-(diff(cuts_tau)*exp(part1+cuts0[1:(nb_index+1)]*alpha[1:(nb_index+1)]))[aInf]
    }
  }
  result<-list(rmst=sum(result),tau=tau)
  class(result)<-"RmstIC"
  return(result)
}
#' @export
print.RmstIC <- function(x, ...)
{
  cat(paste(" restricted mean with upper limit = ",x$tau,"\n",sep=""))
  print(round(x$rmst[1],4))#x$rmst[1]
}
