#' Computation of the maximum likelihood estimator in the PCH model for interval-censored data
#'
#' Given a sequence of time, a set of cuts and initial log-hazard values, compute
#' the maximum likelihood estimator (MLE) of the hazard function from a piecewise-constant hazard model for interval-censored
#' data. The estimator is based on the EM algorithm.
#' @param Left the time sequence for the left-time endpoint. Must be non-negative.
#' @param Right the time sequence for the right-time endpoint. Each component must be greater than the corresponding
#' component of \code{Left}. The value \code{Inf} for right-censored data is allowed.
#' @param cuts the sequence of cuts in the pch model. Default value equals \code{NULL} which corresponds to the exponential model.
#' @param a the initial value of the log-hazard function between each cut. Must be of length equal to \code{length(cuts)+1}.
#' The default value is \eqn{log(0.5)} for all segments.
#' @param maxiter the total maximum number of iterations in the EM algorithm. Default is 1000.
#' @param tol the tolerance value in the EM algorithm. The algorithm stops when the relative error
#' between new and previous log-hazard update value is less than \code{tol}. Default is 1e-12.
#' @param verbose if \code{TRUE}, at each iteration step the current value of the estimator is displayed along with the relative error
#' between new and previous log-hazard update value in the EM algorithm.
#' @details Returns the MLE for the piecewise constant hazard model defined in the following way:
#' \deqn{\lambda(t)=\sum 1(c_{k-1}< t\le c_k) exp(a_k),}
#' where the sum is over \eqn{k} ranging from \eqn{1} to \eqn{K}, with \eqn{K} equals to the length of \code{cuts} plus one,
#' \eqn{c_0=0}, \eqn{c_K} equals infinity and the \eqn{c_k}'s are provided by \code{cuts}. The algorithm estimates the \eqn{a_k}'s using the EM
#' algorithm. The \eqn{a_k}'s are initialized by the value given in \code{a}. The EM algorithm is not sensitive to the initialization
#' and we recommend to simply use the default value which is equal to \eqn{log(0.5)} for each \eqn{a_k}.
#'
#' The algorithm deals with interval-censored data supplied by the \code{Left} and \code{Right} vectors. It is possible
#' to include right-censored data by specifying \code{Right=Inf} and to have exact observations by specifying
#' the same value to \code{Left} and \code{Right}.
#'
#' Two regularity conditions arise from maximum likelihood theory. The first one states that the probability that
#' the observed intervals (for which the right endpoint is not infinite) intersect a cut should be positive. If that is not the case, the
#' algorithm will still be stable, returning the hazard value \eqn{0} for the corresponding cut. The second condition states that the probability
#' that a left endpoint is greater than the left value of a cut should be positive. If there are no observations verifying this
#' condition then the algorithm will diverge and the hazard value for this cut will increase at each iteration of the EM algorithm. This will make the
#' algorithm return an error. More generally, one should carefully check the convergence of the algorithm by setting \code{verbose} to \code{TRUE}.
#'
#' The function \code{mle.ic} returns an object from the \code{survIC} class. For such class, there are a \code{plot} and \code{lines} options.
#' When applying \code{plot} (or \code{lines}) to a \code{survIC} class object, one can choose the option \code{surv=FALSE} (the default)
#' or \code{surv=TRUE}. In the former case, the hazard function will be plotted and in the latter case the survival function will be displayed.
#' One can set \code{xlim} to specific values in the \code{plot} arguments (see the Examples section for more information).
#' @return
#' \tabular{lll}{
#' \code{a} \tab  \code{ } \tab the estimated value of the log-hazard (the \eqn{a_k}'s) \cr
#' \tab \tab \cr
#' \code{lambda} \tab \code{ }  \tab the estimated value of the hazard (the \eqn{exp(a_k)}'s)\cr
#' \tab \tab \cr
#' \code{cuts} \tab  \code{ } \tab the cuts used in the PCH model as initially supplied in the \code{mleIC} function\cr
#' \tab \tab \cr
#' \code{itermax} \tab \code{ }  \tab the number of iterations used in the EM algorithm to reach convergence\cr
#' \tab \tab \cr
#' \code{test} \tab  \code{ } \tab the relative error that has been reached at convergence of the algorithm\cr
#' }
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
#' result=mleIC(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' result$lambda #the estimation
#' alpha #true value
#'
#' plot(result)
#' lines(c(0,cuts,100),c(alpha,alpha[5]),type="s",col="red") #true hazard
#' seqtime=seq(0,200,length.out=1000)
#' plot(result,xlim=c(0,200),surv=TRUE)
#' lines(seqtime,exp(-pchcumhaz(seqtime,cuts,alpha)),type="l",col="red") #true survival function
#' resultbis=mleIC(Left,Right,cuts=NULL,a=rep(log(0.5),1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' lines(resultbis,xlim=c(0,200),col="blue",surv=TRUE)
#'
#' ##Example with exact values
#' Right[1:100]<-TrueTime[1:100]
#' Left[1:100]<-TrueTime[1:100]
#' result=mleIC(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
#' maxiter=1000,tol=1e-12,verbose=TRUE)
#' alpha #true value

#' @export
mleIC=function(Left,Right,cuts=NULL,a=rep(log(0.5),length(cuts)+1),maxiter=1000,tol=1e-12,verbose=FALSE){
  if (length(Left)!=length(Right)) stop("'Left' and 'Right' must have the same length")
  if (all(diff(cuts)>0)==FALSE) {
    stop("cuts must be an increasing sequence with no ex-aequo!")}
  n=length(Left)

  #code to include exact observations
  idEx<-(Left==Right)
  nEx=sum(idEx)
  noEx=n-nEx
  RC<-Right[idEx==FALSE]==Inf
  #RC<-Right==Inf

  test=rep(NA,maxiter)

  cutsInf=c(0,cuts,Inf)
  k=length(cuts)+1
  if (length(a)!=k) stop("'a' must be one length longer than 'cuts'")

  #code to include exact observations
  Idx=lapply(1:k,function(x){Left[idEx==FALSE]<cutsInf[x+1] & Right[idEx==FALSE]>cutsInf[x]})
  Lcl=lapply(1:k,function(x){pmax(cutsInf[x],Left[idEx==FALSE])})
  Rcl=lapply(1:k,function(x){pmin(cutsInf[x+1],Right[idEx==FALSE])})

  #Idx=lapply(1:k,function(x){Left<cutsInf[x+1] & Right>cutsInf[x]})
  #Lcl=lapply(1:k,function(x){pmax(cutsInf[x],Left)})
  #Rcl=lapply(1:k,function(x){pmin(cutsInf[x+1],Right)})

  #code to include exact observations
  if(nEx==0){
    Aex<-Bex<-J<-rep(0,k)} else {
      exact=Left[idEx]
      weights=rep(1,nEx)
      J=cut(exact,breaks=cutsInf,labels=1:k,right=FALSE)
      cuts0=c(0,cuts)
      Aex=sapply(split(weights,J),sum)
      riskFun=function(x){
        vect=c(diff(cuts0),0)
        vect[J[x]]<-(exact[x]-cuts0[J[x]])
        if (as.numeric(J[x])<k){vect[(as.numeric(J[x])+1):k]<-0};return(vect)}
      Bex=t(sapply(1:nEx,riskFun))
      if (k==1){
        Bex=sum(Bex)
      } else {
        Bex=apply(Bex,2,sum)
      }
      }
  for (iter in 1:maxiter)
  {
    old.a=a

    stat=exhaust.ic_No_Cov(n,noEx,Left[idEx==FALSE],Right[idEx==FALSE],cuts,a,Idx,RC,Lcl,Rcl)#exhaust.ic(Left,Right,cuts,a,p)
    APi=stat$APi;BPi=stat$BPi;APiall=stat$APiall;BPiall=stat$BPiall #N=stat$N;
    result=estime_No_Cov(APi,BPi,Aex,Bex,cuts,n)
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
#noEx represents the number of individuals where Left and Right are not equal
exhaust.ic_No_Cov=function(n,noEx,Left,Right,cuts,a,Idx,RC,Lcl,Rcl){
  k=length(cuts)+1
  lambda=exp(a)
  cutsInf=c(0,cuts,Inf)

  cumL=pch.cumhaz(Left,cuts=cuts,a=a)
  cumR=pch.cumhaz(Right,cuts=cuts,a=a)
  cumcl=c(0,pch.cumhaz(cuts,cuts=cuts,a=a))

  if(noEx==0){
    APiall=BPiall=matrix(0,n,k);APi=BPi=rep(0,k)} else {
  logA<-logB<-matrix(-Inf,noEx,k)
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
  APiall=exp(logA)
  BPiall=exp(logB)
  APi=apply(APiall,2,sum)
  BPi=apply(BPiall,2,sum)
  }
  return(list(APi=APi,BPi=BPi,APiall=APiall,BPiall=BPiall))
}

#Return mle estimates from the exhaustive statistics
estime_No_Cov=function(APi,BPi,Aex,Bex,cuts,n){
  #idx=((APi)!=0))
  #code to include exact observations
  idx=((APi+Aex)!=0)
  cumA=c(rev(cumsum(rev(APi[-1]))),0)
  a=rep(-Inf,length(APi))
  a[idx]=log(APi[idx]+Aex[idx])-log(cumA[idx]*c(diff(c(0,cuts)),0)[idx]+BPi[idx]+Bex[idx])
  return(list(a=a))
}

#Compute the objective function Qtheta given thetaold in the EM algorithm
#Not needed
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

