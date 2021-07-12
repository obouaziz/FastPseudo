#' Computation of the maximum likelihood estimator in the PCH model for interval-censored data
#'
#' Given a sequence of time, a set of cuts and initial log-hazard values, compute
#' the maximum likelihood estimator (MLE) of the hazard function from a piecewise-constant hazard model.
#' The estimator is based on the EM algorithm.
#' @param Left the time sequence for the left-time endpoint. Must be non-negative
#' @param Right the time sequence for right-time endpoint. The value \code{Inf} for right-censored data is allowed.
#' @param cuts the sequence of cuts in the pch model. Default is \code{NULL} which corresponds to the exponential model.
#' @param a the initial value of the log-hazard function between each cut. Should be of length equal to \code{length(cuts)+1}.
#' @param maxiter the total number of iterations in the EM algorithm. Default is 1000.
#' @param tol the tolerance value in the EM algorithm. The algorithm stops when the relative error
#' between new and previous log-hazard estimate value is less than \code{tol}. Default is 1e-5.
#' @param verbose if \code{TRUE}, at each iteration step the current value of the estimator is displayed along with the relative error
#' between new and previous log-hazard estimate value.
#' @return return the MLE
#' @export
#' @examples
#' n=4000
#' cuts=c(20,40,50,70)
#' alpha=c(0,0.05,0.1,0.2,0.4)/10
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
#' result=mle.ic(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' result$lambda #the estimation
#' alpha #true value
#' plot(result)
#' lines(c(0,cuts,100),c(alpha,alpha[5]),type="s",col="red") #true hazard
#' seqtime=seq(0,200,length.out=1000)
#' plot(result,xlim=c(0,200),surv=TRUE)
#' lines(seqtime,exp(-pchcumhaz(seqtime,cuts,alpha)),type="l",col="red") #true survival function
#' resultbis=mle.ic(Left,Right,cuts=NULL,a=rep(log(0.5),length(cuts)+1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' lines(resultbis,xlim=c(0,200),col="blue",surv=TRUE)

#' @export
mle.ic=function(Left,Right,cuts=NULL,a,maxiter=1000,tol=1e-5,verbose=FALSE){
  n=length(Left)
  test=Qerr=rep(NA,maxiter)
  RC<-Right==Inf
  cutsInf=c(0,cuts,Inf)
  k=length(cuts)+1

  Idx=lapply(1:k,function(x){Left<cutsInf[x+1] & Right>cutsInf[x]})
  Lcl=lapply(1:k,function(x){pmax(cutsInf[x],Left)})
  Rcl=lapply(1:k,function(x){pmin(cutsInf[x+1],Right)})

  for (iter in 1:maxiter)
  {
    old.a=a

    stat=exhaust.ic_No_Cov(n,Left,Right,cuts,a,Idx,RC,Lcl,Rcl)#exhaust.ic(Left,Right,cuts,a,p)
    N=stat$N;APi=stat$APi;BPi=stat$BPi;APiall=stat$APiall;BPiall=stat$BPiall
    result=estime_No_Cov(APi,BPi,cuts,n)
    a=result$a
    a_noInf<-a!=-Inf
    test[iter]=max(abs(a[a_noInf]-old.a[a_noInf])/abs(a[a_noInf]))
    if (verbose==TRUE){
      cat("lambda=",exp(a),"abs.err=",test[iter],"iter=",iter,"\n",sep=" ")}#test[iter]
    if (test[iter]<tol) break;
  }
  hazard=NULL
  if (k>1) hazard=stats::stepfun(cuts,exp(a))
  out<-list(a=a,cuts=cuts,lambda=exp(a),hazard=hazard,APi=APi,BPi=BPi,APiall=APiall,BPiall=BPiall,itermax=iter,test=test[1:iter])
  class(out)<-"survIC"
  out
}

#' @export
plot.survIC<-function(x,...,xlab="Time",ylab="Survival",main="Survival function",xlim=NULL,surv=FALSE)
{
  if (surv==FALSE){
    if (is.null(x$cuts)) {stop("There is no plot for the exponential hazard model")} else {
      stats::plot.stepfun(x$hazard,do.points=FALSE,xlab="Time",ylab="Hazard",main="Hazard function",...)}
      } else {
        if (is.null(x$cuts)) {
          if (is.null(xlim)) {stop("Please provide an xlim argument")} else {
      seqtime=seq(xlim[1],xlim[2],length.out=1000)
      graphics::plot.default(seqtime,exp(-pch.cumhaz(seqtime,x$cuts,x$a)),type="l",xlab=xlab,ylab=ylab,main=main,...)
          }
        } else {
          if (is.null(xlim)) {xlim=c(0,max(x$cuts)+max(diff(c(0,x$cuts))))}
          seqtime=seq(xlim[1],xlim[2],length.out=1000)
          graphics::plot.default(seqtime,exp(-pch.cumhaz(seqtime,x$cuts,x$a)),type="l",xlab=xlab,ylab=ylab,main=main,...)
          }
        }
}

#' @export
lines.survIC=function(x,...,xlim=NULL,surv=FALSE)
{
  if (surv==FALSE){
    if (is.null(x$cuts)) {graphics::abline(h=x$lambda,...)} else {
      if (is.null(xlim)) {xlim=c(0,max(x$cuts)+max(diff(c(0,x$cuts))))}
      graphics::lines.default(c(0,x$cuts,xlim[2]),c(x$lambda,x$lambda[length(x$cuts)+1]),type="s",...)
  }} else {
    if (is.null(x$cuts)) {
      if (is.null(xlim)) {stop("Please provide an xlim argument")} else {
        seqtime=seq(xlim[1],xlim[2],length.out=1000)
    graphics::lines.default(seqtime,exp(-pch.cumhaz(seqtime,x$cuts,x$a)),type="l",...)
      }
    } else {
      if (is.null(xlim)) {xlim=c(0,max(x$cuts)+max(diff(c(0,x$cuts))))}
      seqtime=seq(xlim[1],xlim[2],length.out=1000)
      graphics::lines.default(seqtime,exp(-pch.cumhaz(seqtime,x$cuts,x$a)),type="l",...)
    }
  }
}

# given cuts and a=log(hazard), return the cumulative hazard function evaluated at the vector time
#same function as pchcumhaz except that it is on the log scale
pch.cumhaz=function(Time,cuts,a) {
  if (length(cuts)==0) return(Time*exp(a[1]))
  k=length(a)
  J=cut(Time,breaks=c(-1,cuts,Inf),labels=1:k,right=TRUE)
  cuts0=c(0,cuts)
  return(exp(a)[J]*(Time-cuts0[J])+c(0,cumsum(exp(a)*c(diff(cuts0),0))[-k])[J])
}


#Return the exhaustive statistics for one loop inside the EM algorithm
exhaust.ic_No_Cov=function(n,Left,Right,cuts,a,Idx,RC,Lcl,Rcl){
  k=length(cuts)+1
  lambda=exp(a)
  cutsInf=c(0,cuts,Inf)

  cumL=pch.cumhaz(Left,cuts=cuts,a=a)
  cumR=pch.cumhaz(Right,cuts=cuts,a=a)
  cumcl=c(0,pch.cumhaz(cuts,cuts=cuts,a=a))

  logA<-logB<-matrix(-Inf,n,k)
  for (l in 1:k)
  { #this condition is useful only if lambda[1]==0 and there exists
    #L_i, R_i s.t L_i<c_1 and R_i<c_1, this gives S(L_i)-S(R_i)=0
    #which is the denominator in A_{k,i},B_{k,i}
    if (lambda[l]!=0){
      id=Idx[[l]];L=Lcl[[l]][id];R=Rcl[[l]][id]
      expterm=exp(-lambda[l]*(R-L))
      lastpart=lambda[l]*cutsInf[l]-cumcl[l]+cumL[id]-log1p(-exp(cumL[id]-cumR[id]))
      logA[id,l]=-lambda[l]*L+log1p(-expterm)+lastpart
      #}
      Bterm=rep(0,sum(id));newid=(R!=Inf)
      if (sum(newid)!=0){
        ratio=(L[newid]-R[newid])/(exp(-a[l])+R[newid]-cutsInf[l])
        Bterm[newid]=log1p(-exp(log(expterm[newid])-log1p(ratio)))}
      logB[id,l]=log(1/lambda[l]+L-cutsInf[l])-lambda[l]*L+Bterm+lastpart
    }
  }
  return(list(APi=apply(exp(logA),2,sum),BPi=apply(exp(logB),2,sum),APiall=exp(logA),BPiall=exp(logB)))
}

#Return mle estimates from the exhaustive statistics
estime_No_Cov=function(APi,BPi,cuts,n){
  idx=(APi!=0)
  cumA=c(rev(cumsum(rev(APi[-1]))),0)
  a=rep(-Inf,length(APi))
  a[idx]=log(APi[idx])-log(cumA[idx]*c(diff(c(0,cuts)),0)[idx]+BPi[idx])
  return(list(a=a))
}

#Compute the objective function Qtheta given thetaold in the EM algorithm
#Probably not needed
qtheta.ic_No_Cov=function(APi,BPi,cuts,n,a,pen=0,w=0) {
  #if (is.null(weights)) weights=rep(1,n)
  k=length(cuts)+1
  lambda=exp(a)
  bracket=a-c(0,cumsum(c(diff(c(0,cuts[-k])))*c(lambda[-k])))
  term2=sum(bracket*APi-lambda*BPi)
  result=term2-pen*sum(w*(diff(a))^2)/2
  return(result)
}

#' #'
#' plot.survIC=function(fit,xlim=c(0,max(fit$cuts)+max(diff(c(0,fit$cuts)))),xlab="Time",ylab="Survival",main="Survival function",lwd=1,...,surv=FALSE)
#' {
#'   if (surv==FALSE){
#'     stats::plot.stepfun(fit$hazard,do.points=FALSE,xlab="Time",ylab="Hazard",main="Hazard function",lwd=lwd,...)
#'   } else {
#'   seqtime=seq(xlim[1],xlim[2],length.out=1000)
#'   graphics::plot.default(seqtime,exp(-pch.cumhaz(seqtime,fit$cuts,fit$a)),type="l",xlab=xlab,ylab=ylab,main=main,lwd=lwd,...)
#'   }
#' }

#' #'
#' lines.survIC=function(fit,xlim=c(0,max(fit$cuts)+max(diff(c(0,fit$cuts)))),lwd=1,...,surv=FALSE)
#' {
#'   if (surv==FALSE){
#'     stats::lines.stepfun(fit$hazard,do.points=FALSE,lwd=lwd,...)
#'   } else {
#'     seqtime=seq(xlim[1],xlim[2],length.out=1000)
#'     graphics::lines.default(seqtime,exp(-pch.cumhaz(seqtime,fit$cuts,fit$a)),type="l",lwd=lwd,...)
#'   }
#' }

