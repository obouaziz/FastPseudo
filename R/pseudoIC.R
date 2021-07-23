#' Fast pseudo values for the survival function or the RMST for interval-censored data
#'
#'  Compute pseudo values for the survival function or the Restricted Mean Survival Time (RMST) based on the piecewise-constant hazard (PCH) model.
#'  A fast approximation is implemented that provides asymptotic pseudo-values based on the score information
#'  and Hessian matrix instead of computing the jackknife method.
#' @param fit an object from the class \code{survIC}. This object is obtained by implementing the maximum likelihood
#' estimator in the PCH model for interval-censored data through the function \code{mleIC}.
#' @param Left the time sequence for the left-time endpoint. Must be non-negative.
#' @param Right the time sequence for the right-time endpoint. Each component must be greater than the corresponding
#' component of \code{Left}. The value \code{Inf} for right-censored data is allowed.
#' @param tseq is used only for the implementation of the pseudo-values of the survival function
#' (in which case tau=\code{NULL}). The survival function is evaluated at the times
#' given in \code{tseq} and its pseudo-values are computed.
#' @param tau is equals \code{NULL} then the program implements the pseudo-values
#' for the survival function. Otherwise, tau must be a positive time point in which
#' case the RMST with right endpoint tau is implemented.
#' @details A fast approximation formula is used to compute pseudo-values for interval-censored data based on the pch model. For more information on
#' the pch model, see the help for the function \code{mleIC}.  The pseudo-values are computed
#' for the survival function if \code{tau} is equal to \code{NULL}. In that case, the function returns an approximation of
#'
#' \deqn{n\hat S(t)-(n-1)\hat S^{(-i)}(t),}
#'
#' where \eqn{\hat S} is the survival estimator from the pch model and \eqn{\hat S^{(-i)}} is the same estimator computed when the \eqn{i}th
#' observation has been removed. Those pseudo-values are computed for all \eqn{i=1,...,n} and for all time points \eqn{t} in \code{tseq}. They are computed for the RMST if \code{tau} is specified. In that case, the function returns an approximation of
#'
#' \deqn{n\int_0^{\tau}\hat S(t)dt-(n-1)\int_0^{\tau}\hat S^{(-i)}(t)dt.}
#'
#' A model based on Generalised Estimating Equation can be further specified through the \code{geese} function
#' in the \code{geepack} package (see Examples below).
#'
#' Regularity conditions are imposed from maximum likelihood theory. Those conditions may not be
#' verified if the cuts are not adequately chosen in the \code{fit} object and  this will be problematic for the computation of the pseudo-values.
#' In particular, the average of the pseudo-values might not return the value
#' of the initial estimator, indicating that the pseudo-values are incorrect. We therefore recommend to always carefully check
#' if the model in \code{fit} implemented from the \code{mleIC} function has properly converged. This can easily be
#' verified by using the option \code{verbose=TRUE} in the \code{mleIC} function. If the \code{mleIC} function returned an error or
#' if the estimate obtained from the \code{mleIC} function contains some \eqn{\alpha_k}'s with the value \eqn{0}, then the
#' \code{pseudoIC} function will return an error. We then recommend
#' to simply change the values of the cuts. See \code{mleIC} for more details about
#' those regularity conditions.
#'
#' @seealso \code{\link{mleIC}}, \code{\link{RmstIC}}, \code{\link{rsurv}}, \code{\link{pseudoKM}}.
#'
#' @return if \code{tau}=\code{NULL}, returns the pseudo-values of the survival function
#' evaluated at \code{tseq}. The output is a matrix with number of rows equal to
#' the sample size and number of columns equal to the length of \code{tseq}.
#'
#' if \code{tau} contains a single positive value then the function returns the
#' RMST with right endpoint equal to \code{tau}. In case \code{tseq} and \code{tau}
#' are both provided, the RMST is implemented and the value of \code{tseq} is ignored.
#'
#' @examples
#' #Simulations without covariates
#' n=4000
#' cuts=c(20,40,50,70)
#' alpha=c(0.05,0.05,0.1,0.2,0.4)/10
#' TrueTime=rsurv(n,cuts,alpha) #generate true data from the pch model
#' ##Simulation of interval-censored data
#' Right<-rep(Inf,n)
#' nb.visit=20
#' visTime=0;visit=matrix(0,n,nb.visit+1)
#' visit=cbind(visit,rep(Inf,n))
#' visit[,2]=visit[,1]+stats::runif(n,0,10)#runif(n,0,5)
#' schedule=4
#' for (i in 3:(nb.visit+1))
#' {
#'   visit[,i]=visit[,i-1]+stats::runif(n,0,schedule*2)
#' }
#' Left<-visit[,(nb.visit+1)]
#' J=sapply(1:(n),function(i)cut(TrueTime[i],breaks=c(visit[1:(n),][i,]),
#' labels=1:(nb.visit+1),right=FALSE)) #sum(is.na(J)) check!
#' Left[1:(n)]=sapply(1:(n),function(i)visit[1:(n),][i,J[i]])
#' Right[1:(n)]=sapply(1:(n),function(i)visit[1:(n),][i,as.numeric(J[i])+1])
#' #View(data.frame(Left,Right,TrueTime)) #To see the generated data
#' sum(Right[1:(n)]==Inf)/n #percentage of right-censored data
#'
#' result=mleIC(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' result$lambda #the estimation
#' alpha #true value of hazard
#' #Plot the true survival function and compare with pseudo-values
#' plot(result,surv=TRUE)
#' tseq=seq(30,90,length.out=10)
#' pseudo_val=pseudoIC(result,Left,Right,tseq=tseq,tau=NULL)
#' lines(tseq,apply(pseudo_val,2,mean),type="p",col="red")
#'
#' #Estimated value of RMST for tau=30
#' RmstIC(cuts,result$lambda,tau=30)
#' pseudoval=pseudoIC(result,Left,Right,tau=30) #pseudo-values
#' mean(pseudoval)
#' #Estimated value of RMST for tau=45
#' RmstIC(cuts,result$lambda,tau=45)
#' pseudoval=pseudoIC(result,Left,Right,tau=45) #pseudo-values
#' mean(pseudoval) #is close to the estimated RMST
#' #Estimated value of RMST for tau=50
#' RmstIC(cuts,result$lambda,tau=50)
#' pseudoval=pseudoIC(result,Left,Right,tau=50) #pseudo-values
#' mean(pseudoval) #is close to the estimated RMST
#'
#' #Simulations with covariates
#' n=1000
#' tau=3000
#' X=runif(n,0,2)
#' epsi=rnorm(n,0,1)
#' theta0=6;theta1=4
#' TrueTime=theta0+theta1*X+epsi #min(TrueTime) must be positive!
#' ##Simulation of interval-censored data
#' Right<-rep(Inf,n)
#' nb.visit=5#nb.visit=10
#' visTime=0;visit=matrix(0,n,nb.visit+1)
#' visit=cbind(visit,rep(Inf,n))
#' visit[,2]=visit[,1]+runif(n,0,10)#runif(n,0,5)
#' schedule=2
#' for (i in 3:(nb.visit+1))
#' {
#'   visit[,i]=visit[,i-1]+runif(n,0,schedule*2)
#' }
#' Left<-visit[,(nb.visit+1)]
#' J=sapply(1:(n),function(i)cut(TrueTime[i],breaks=c(visit[1:(n),][i,]),
#' labels=1:(nb.visit+1),right=FALSE)) #sum(is.na(J)) check!
#' Left[1:(n)]=sapply(1:(n),function(i)visit[1:(n),][i,J[i]])
#' Right[1:(n)]=sapply(1:(n),function(i)visit[1:(n),][i,as.numeric(J[i])+1])
#'
#' cuts=c(6,8,10,12,14)
#' result=mleIC(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' pseudoval=pseudoIC(result,Left,Right,tau=tau)
#' require(geepack)
#' data_pseudo<-data.frame(Y=c(pseudoval),X=X,id=rep(1:n))
#' resultEst=geese(Y~X,id=id,jack=TRUE,family="gaussian", mean.link="identity",
#' corstr="independence",scale.fix=TRUE,data=data_pseudo)
#' summary(resultEst)





#' @export
pseudoIC<-function(fit,Left,Right,tseq=NULL,tau=NULL){
  if (length(Left)!=length(Right)) stop("'Left' and 'Right' must have the same length")
  cuts=fit$cuts
  APi=fit$APi
  BPi=fit$BPi
  a=fit$a
  if (min(a)==-Inf) warning("Some values of the hazard are equal to 0.
                            The pseudo-values will not be valid, please change
                            the value of the cuts in the fit object")
  n=length(Left)
  k=length(cuts)+1
  O_lIC=fit$APi
  cumA=c(rev(cumsum(rev(APi[-1]))),0)
  R_lIC=cumA*c(diff(c(0,cuts)),0)+BPi
  NoRC=Right!=Inf
  num=t(sapply(Right[NoRC],function(x) diffcut(x,cuts))-sapply(Left[NoRC],function(x) diffcut(x,cuts)))
  if (k==1) {num=t(num)}
  x=pch.cumhaz(Left[NoRC],cuts,a)-pch.cumhaz(Right[NoRC],cuts,a)

  term2=log1p(-exp(x))
  term2=x-2*term2
  Hess=t(num)%*%(num*matrix(rep(exp(term2),k),nrow=sum(NoRC)))/n

  num2=t(sapply(Right,function(x) diffcut(x,cuts))-sapply(Left,function(x) diffcut(x,cuts)))
  if (k==1) {num2=t(num2)}
  num2[Right==Inf,]<-rep(0,k)
  x2=pch.cumhaz(Left,cuts,a)-pch.cumhaz(Right,cuts,a)
  term3=x2-log1p(-exp(x2))
  term3[Right==Inf]<-0
  Score_part1=-t(sapply(Left,function(x) diffcut(x,cuts)))
  if (k==1){Score_part1=t(Score_part1)}
  Score2=Score_part1+num2*matrix(rep(exp(term3),k),nrow=n)#sum(NoRC)

  OR_frac=t(solve(Hess,t(Score2)))
  if (is.null(tau)) #Compute the pseudo-value for the pch survival function
  {
    if (is.null(tseq)) stop("Either a value for 'tau' or for 'tseq' must be provided")
    if (length(tseq)==1){
      surv=exp(-pch.cumhaz(tseq,cuts = cuts,a = a))
      gradLambda=as.matrix(diffcut(tseq,cuts))
    } else {
      surv=matrix(rep(exp(-pch.cumhaz(tseq,cuts = cuts,a = a)),n),byrow=TRUE,nrow=n)
      gradLambda= sapply(tseq,function(x){diffcut(x,cuts)})
    }
    pseudoval=surv-surv*OR_frac%*%gradLambda
  } else { #compute the pseudo-value for the RMST in the pch model
    if (is.null(tseq)==FALSE) warning("The value of tseq is ignored and the RMST is returned")
  cuts0=c(0,cuts)
  if (k==1)
  {
    pseudo_Est=resM(cuts,a,tau)$surv_int-OR_frac*resM2(cuts,a,tau)$surv_int
  } else {
    if(length(cuts)==1){
      B=cbind(rep(0,n),OR_frac[,-k]*diff(cuts0))
    } else {
      B=cbind(rep(0,n),t(apply(OR_frac[,-k]*matrix(rep(diff(cuts0),n),byrow=TRUE,nrow=n),1,cumsum)))
    }
    #cumsum only on lines such that the sum is on k
    B=B-matrix(rep(cuts0,n),byrow=TRUE,nrow=n)*OR_frac #add the last piece of B
    index=c(tau>cuts) #watch out that tau>0!c(TRUE,tau>cuts)
    nb_index=sum(index)
    B=B[,1:(nb_index+1)]
    part1=B*matrix(rep(resM(cuts,a,tau)$surv_int_cut,n),byrow=TRUE,nrow=n)
    part2=OR_frac[,1:(nb_index+1)]*matrix(rep(resM2(cuts,a,tau)$surv_int_cut,n),byrow=TRUE,nrow=n)
    pseudoval=resM(cuts,a,tau)$surv_int-apply(part1+part2,1,sum)
  }
  }
  return(pseudoval)
}


#Compute c_k\wedge\tau-c_{k-1} where obsval is tau.
#obsval must be a unique value (cannot be a vector)!
diffcut=function(obsval,cuts){
  k=length(cuts)+1
  J=cut(obsval,breaks=c(-1,cuts,Inf),labels=1:k,right=TRUE)
  cuts0=c(0,cuts)
  if (as.numeric(J)==1){return(c(obsval,rep(0,k-1)))} else{
    return(c(c(diff(cuts0))[1:(as.numeric(J)-1)],obsval-cuts0[J],rep(0,k-as.numeric(J))))
  }
}

#Compute the RSMT for a fixed value of tau
resM=function(cuts=NULL,a,tau){
  alpha=exp(a)
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
  return(list(surv_int_cut=result,surv_int=sum(result)))
}

#Compute the integral of t*S(t) from 0 to tau
resM2=function(cuts=NULL,a,tau){
  alpha=exp(a)
  if (is.null(cuts))
  {
    result=exp(log1p(-exp(-alpha*tau)*(1+alpha*tau))-2*a)
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
    part2=log1p(-exp(-alpha[1:(nb_index+1)]*diff(cuts_tau))*(1+alpha[1:(nb_index+1)]*cuts_tau[-1])+alpha[1:(nb_index+1)]*cuts_tau[-(nb_index+2)])-2*a[1:(nb_index+1)]
    result=(exp(part1+part2))
    if (sum(a[1:(nb_index+1)]==-Inf)>0) {
      aInf<-(a[1:(nb_index+1)]==-Inf)
      result[aInf]<-((((cuts_tau[-1])^2-(cuts0[1:(nb_index+1)])^2))*exp(part1+cuts0[1:(nb_index+1)]*alpha[1:(nb_index+1)])/2)[aInf]###to check!
    }
  }
  return(list(surv_int_cut=result,surv_int=sum(result)))
}

#Estimated value of RMST for tau=45
# 20+exp(alpha[2]*20)*(exp(-alpha[2]*20)-exp(-alpha[2]*30))/alpha[2]
# #true value of RMST for tau=30
# 20+exp(result$lambda[2]*20)*(exp(-result$lambda[2]*20)-
# exp(-result$lambda[2]*30))/result$lambda[2]
# #estimated value of RMST for tau=30
# pseudoval=pseudoIC(result,Left,Right,tau=30)
# mean(pseudoval) #is close to the estimated RMST
# #true value of RMST for tau=45
# 20+exp(result$lambda[2]*20)*(exp(-result$lambda[2]*20)-
# exp(-result$lambda[2]*40))/result$lambda[2]+
# exp(-20*result$lambda[2]+40*result$lambda[3])*(exp(-result$lambda[3]*40)-
# exp(-result$lambda[3]*45))/result$lambda[3]
# #estimated value of RMST for tau=45
# 20+exp(alpha[2]*20)*(exp(-alpha[2]*20)-exp(-alpha[2]*40))/alpha[2]+
# exp(-20*alpha[2]+40*alpha[3])*(exp(-alpha[3]*40)-exp(-alpha[3]*45))/alpha[3]
# pseudoval=pseudoIC(result,Left,Right,tau=45)
# mean(pseudoval) #is close to the estimated RMST
# #' #true value of RMST for tau=50
# 20+exp(result$lambda[2]*20)*(exp(-result$lambda[2]*20)-
# exp(-result$lambda[2]*40))/result$lambda[2]+
# exp(-20*result$lambda[2]+40*result$lambda[3])*(exp(-result$lambda[3]*40)-
# exp(-result$lambda[3]*50))/result$lambda[3]
# pseudoval=pseudoIC(result,Left,Right,tau=50)
# mean(pseudoval) #is close to the estimated RMST




