#' Fast pseudo values for the survival function and the RMST for right-censored data
#'
#'  Compute pseudo values for the Kaplan-Meier estimator or for the Restricted Mean Survival Time (RMST) based on the Kaplan-Meier estimator.
#'  We use a fast approximation of the jackknife method.
#'
#' @param Time a continuous time variable.
#' @param status an indicator for censoring for the corresponding Time variable.
#' @param tau if equal to NULL, the pseudo-values for the Kaplan-Meier estimator are computed. If a value of tau is
#' specified, the pseudo-values for the RMST with time endpoint equal to tau is computed.
#' @details A fast approximation formula is used to compute pseudo-values for right-censored data. Those pseudo-values are computed
#' for the survival function if \code{tau} is equal to \code{NULL}. In that case, the function returns an approximation of
#'
#' \deqn{n\hat S(t)-(n-1)\hat S^{(-i)}(t),}
#'
#' where \eqn{\hat S} is the Kaplan-Meier estimator and \eqn{\hat S^{(-i)}} is the Kaplan-Meier estimator computed when the \eqn{i}th
#' observation has been removed. Those pseudo-values are computed for all \eqn{i=1,...,n} and for all time points \eqn{t} corresponding
#' to all observed times. Those pseudo-values are computed for the RMST if \code{tau} is specified. In that case, the function returns an approximation of
#'
#' \deqn{n\int_0^{\tau}\hat S(t)dt-(n-1)\int_0^{\tau}\hat S^{(-i)}(t)dt.}
#'
#' A model based on Generalised Estimating Equation can be further specified through the \code{geese} function
#' in the \code{geepack} package (see Examples below).
#' @return For the survival function (when tau equals NULL) it returns the pseudo values in a matrix form where the individuals are in column and the times
#' are in the rows. For the RMST (when tau is specified) it returns the pseudo values in vector form.
#' @export
#' @examples
#' #Illustration of the pseudo-values for the Kaplan-Meier estimator on a simple simulated set
#' require(survival)
#' n=100
#' cpar=0.1
#' set.seed(28)
#' TrueTime=rweibull(n,shape=0.5,scale=2)
#' Cens=rexp(n,cpar)
#' Tobs=pmin(TrueTime,Cens) #observed times
#' status=TrueTime<=Cens #mean(status) #21% of censoring on average
#'
#' pseudo_val=pseudoKM(Tobs,status,tau=NULL)
#' VonM=pseudo_val$pseudoval
#' tseq=pseudo_val$tseq
#'
#' par(mfrow=c(2,2))
#' plot(tseq,VonM[,3],type="s",xlab="Time",ylab="Pseudo-val. for obs. 3")
#' abline(v=Tobs[3],lty=2,col="red")
#' plot(tseq,VonM[,48],type="s",xlab="Time",ylab="Pseudo-val. for obs. 48")
#' abline(v=Tobs[48],lty=2,col="red")
#' plot(tseq,VonM[,1],type="s",xlab="Time",ylab="Pseudo-val. for obs. 1")
#' abline(v=Tobs[1],lty=2,col="red")
#' plot(tseq,VonM[,5],type="s",xlab="Time",ylab="Pseudo-val. for obs. 5")
#' abline(v=Tobs[5],lty=2,col="red")
#' par(mfrow=c(1,1))
#'
#' plot(tseq,apply(VonM,1,mean),type="s",xlab="Time",ylab="Mean of pseudo-values")
#' lines(survfit(Surv(Tobs,status)~1),col="red",conf.int=FALSE)
#'
#' #Illustration on simple simulated data for the pseudo-values for the RMST
#' #Simulation in the Exponential model (no covariates)
#'
#' n<-100
#' lambda=0.5
#' TrueTime=rexp(n,rate=lambda)
#' tau=10
#' tau=c(2,4,6,8,10)
#' cpar=0.2
#' Cens=rexp(n,cpar)
#' Tobs=pmin(TrueTime,Cens) #observed times
#' status=TrueTime<=Cens #mean(status) #36% of censoring on average
#'
#' pseudo_val=sapply(tau,function(x){pseudoKM(Tobs,status,tau=x)$pseudoval})
#' apply(pseudo_val,2,mean) #check that mean of pseudo-values return the RMST
#' sapply(tau,function(x){Rmst(Tobs,status,x)$rmst}) #compute the RMST
#' sapply(tau,function(x){(1-exp(-lambda*x))/lambda}) #the truth
#'
#'#Illustration on simple simulated data for the pseudo-values for the RMST
#'#Simulation in a linear model
#' set.seed(28)
#' require(geepack)
#' n<-4000
#' sigma=1
#' tau=400
#' alpha=c(3)
#' X=rnorm(n,2,0.25)
#' epsi=rnorm(n,0,sigma)
#' TrueTime=alpha*X+epsi
#' Cens=rep(Inf,n)

#' Tobs=pmin(TrueTime,Cens)
#' status=TrueTime<=Cens
#'
#' VonM=pseudoKM(Tobs,status,tau)$pseudoval
#' data_VonM<-data.frame(Y=c(VonM),X=X,id=rep(1:n))
#' resultEst=geese(Y~X,id=id,jack=TRUE,family="gaussian", mean.link="identity",
#' corstr="independence",scale.fix=TRUE,data=data_VonM)
#' summary(resultEst)
#'
#'
#' #Illustration on simulated data for the estimation in a Cox model based on
#' #the pseudo-values of the Kaplan-Meier estimator
#' #Compare pseudoKM with jackknife method
#'
#' set.seed(28)
#' require(geepack)
#' require(survival)
#' n=4000
#' shape=3;scale=30;cpar=0.01
#' #Simulate covariates
#' X=runif(n,0,1)
#' #simulate Cox model with Weibull baseline distribution
#' TrueTime=scale*(rexp(n,1)/exp(log(0.5)*X))^(1/shape)
#' #True regression parameter is log(0.5)=-0.6931472
#' Cens=rexp(n,cpar)
#' Tobs=pmin(TrueTime,Cens) #observed times
#' Tsort<-sort(Tobs,index.return=TRUE)
#' Tobs_ord<-Tsort$x
#' status=TrueTime<=Cens #mean(status) #25% of censoring on average
#' status_ord<-status[Tsort$ix]
#' X_ord=X[Tsort$ix]
#' SurvEst<-matrix(NA,n,n)
#' #The jackknife method (takes some time)
#' for (i in 1:n)
#' {
#'   Tjack<-Tobs_ord[-i]
#'   statusjack<-status_ord[-i]
#'   resKMjack=survfit(Surv(Tjack,statusjack)~1)
#'   resKMfunjack=stepfun(c(resKMjack$time),c(1,resKMjack$surv))
#'   SurvEst[,i]<-resKMfunjack(c(Tobs_ord))#[status_ord==1]
#' }
#' #Compute pseudo-values
#' resKM=survfit(Surv(Tobs,status)~1)
#' pseudo<-n*matrix(rep(resKM$surv,n),ncol=n)-(n-1)*SurvEst
#' pseudo=pseudo[n*seq(5,95,by=10)/100,] #we compute the pseudo-values for 10 times
#' tseq=Tobs_ord[n*seq(5,95,by=10)/100]
#' M=length(n*seq(5,95,by=10)/100)
#' data_pseudo<-data.frame(Y=1-c(pseudo),X=rep(X_ord,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
#' result=geese(Y~X+as.factor(Time)-1,id=id,jack=TRUE,
#' family="gaussian", mean.link="cloglog", corstr="independence",scale.fix=TRUE,data=data_pseudo)
#' summary(result)
#' ###########################################
#' #with pseudoKM
#' VonM=pseudoKM(Tobs,status,tau=NULL)$pseudoval
#' pseudo_VM=VonM[n*seq(5,95,by=10)/100,] #we compute the pseudo-values for 10 times
#' tseq=Tobs_ord[n*seq(5,95,by=10)/100]
#' M=length(n*seq(5,95,by=10)/100)
#' data_pseudo1<-data.frame(Y=1-c(pseudo_VM),X=rep(X,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
#' #data_pseudo1<-data.frame(Y=1-c(pseudo_VM),X=rep(X_ord,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
#' result1<-geese(Y~X+as.factor(Time)-1,id=id,jack=TRUE,family="gaussian",
#' mean.link="cloglog", corstr="independence",scale.fix=TRUE,data=data_pseudo1)
#' summary(result1)
#'
#' #Illustration on simulated data for the pseudo-values for the RMST
#'#Simulation in a linear model
#' set.seed(28)
#' n<-10000
#' sigma=3
#' cpar=0.07
#' tau=4
#' alpha=c(5.5,0.25,0.25)
#' X1=rbinom(n,1,0.5)
#' X2=rbinom(n,1,0.5)
#' epsi=runif(n,-sigma,sigma)
#' TrueTime=alpha[1]+alpha[2]*X1+alpha[3]*X2+epsi
#' Cens=rexp(n,cpar)

#' Tobs=pmin(TrueTime,Cens)
#' Tsort<-sort(Tobs,index.return=TRUE)
#' Tobs_ord<-Tsort$x
#' status=TrueTime<=Cens #approximately 35% of censoring
#' status_ord<-status[Tsort$ix]
#' X_ord1=X1[Tsort$ix];X_ord2=X2[Tsort$ix]

#' X00=(X1==0 & X2==0)
#' X01=(X1==0 & X2==1)
#' X10=(X1==1 & X2==0)
#' X11=(X1==1 & X2==1)
#' X_ord00=X00[Tsort$ix];X_ord01=X01[Tsort$ix];X_ord10=X10[Tsort$ix];X_ord11=X11[Tsort$ix]
#' #The true value of the parameters are c(3.812552,0.05705492,0.05730445,0.1044318)
#' #Those values were based on a Monte-Carlo experiment with n=10 000 000, using the true times
#'
#' VonM=pseudoKM(Tobs,status,tau)$pseudoval
#' data_VonM<-data.frame(Y=c(VonM),X2=X01,X3=X10,X4=X11,id=rep(1:n))
#' resultEst=geese(Y~X2+X3+X4,id=id,jack=TRUE,family="gaussian",
#' mean.link="identity", corstr="independence",scale.fix=TRUE,data=data_VonM)
#' summary(resultEst)




#' @export
pseudoKM <- function(Time,status,tau=NULL) {UseMethod("pseudoKM")}
#' @export
pseudoKM.default<-function(Time,status,tau=NULL){
  if (length(Time)!=length(status)) stop("'Time' and 'status' must have the same length")
  Tsort<-sort(Time,index.return=TRUE)
  Tobs_ord<-Tsort$x
  status_ord<-status[Tsort$ix]
  if (is.null(tau))#(rmst==FALSE)
  {
    pseudoval=pseudo_ord(Tobs_ord,status_ord)
    pseudoval=pseudoval[,order(Tsort$ix)] #put the values in initial order
    #tseq=Tobs_ord
    result<-list(pseudoval=pseudoval,tseq=Tobs_ord)

  } else {
  #if (rmst==TRUE)
  {
    pseudoval=pseudoRMST_ord(Tobs_ord,status_ord,tau)
    pseudoval=pseudoval[order(Tsort$ix)] #put the values in initial order
    result<-list(pseudoval=pseudoval)
  }
  }
  #result<-list(pseudoval=pseudoval,tau=tau)
  class(result)<-"pseudoKM"#result
  return(result)#result
}

#' @export
pseudoKM.Surv <- function(Time,status=NULL,tau=NULL)#,alternative="two.sided"
{
  formula=Time
  tau2=tau
  Time<-formula[,1]
  status<-formula[,2]
  Tsort<-sort(Time,index.return=TRUE)
  Tobs_ord<-Tsort$x
  status_ord<-status[Tsort$ix]
  if (is.null(tau2))#(rmst==FALSE)
  {
    pseudoval=pseudo_ord(Tobs_ord,status_ord)
    pseudoval=pseudoval[,order(Tsort$ix)] #put the values in initial order
    #tseq=Tobs_ord
    result<-list(pseudoval=pseudoval,tseq=Tobs_ord)

  } else {
    #if (rmst==TRUE)
    {
      pseudoval=pseudoRMST_ord(Tobs_ord,status_ord,tau2)
      pseudoval=pseudoval[order(Tsort$ix)] #put the values in initial order
      result<-list(pseudoval=pseudoval)
    }
  }
  #result<-list(pseudoval=pseudoval,tau=tau)
  class(result)<-"pseudoKM"#result
  return(result)#result
}

### Helper ####
#Compute the pseudo-values for the survival function for RC data for all time points corresponding to observations
#Works for ordered times only
pseudo_ord<-function(Tobs_ord,status_ord){
  n<-length(Tobs_ord)
  #if (n!=length(status_ord)) stop("'Time' and 'status' must have the same length")
  resKM=survival::survfit(survival::Surv(Tobs_ord,status_ord)~1)
  VonM<-matrix(NA,n,n)
  H<-seq(1,1/n,by=-1/n)
  for (l in 1:n)
  {
    part1=cumsum(status_ord/(((n-seq(1,n,by=1)+1)^2)/n))
    part1=c(part1[1:l],rep(part1[l],n-l))
    VonM[,l]=resKM$surv*(1+part1-status_ord[l]*(1/H[l])*c(rep(0,l-1),rep(1,n-l+1)))
  }
  return(VonM)
}

#Compute the pseudo-values for RMST for RC data with endpoint specified by tau
#Works for ordered times only
pseudoRMST_ord<-function(Tobs_ord,status_ord,tau){
  #if (is.na(tau)) stop("'tau' cannot be equal to NA")
  #if (tau==Inf) stop("'tau' cannot be infinite")
  if (is.na(tau)|(tau==Inf)) stop("a finite value of tau must be specified")
  n<-length(Tobs_ord)
  if (n!=length(status_ord)) stop("'Time' and 'status' must have the same length")
  if (sum(Tobs_ord<=tau)<3) stop("the value of tau is too small")
  VonM<-rep(NA,n)
  resKM=survival::survfit(survival::Surv(Tobs_ord,status_ord)~1)
  #for large sample sizes, sometimes some Tobs are removed in survfit (so in resKM)...
  #then length(resKM$surv) is <n and the codes won't work
  #this is a fix
  if (length(resKM$surv)!=n){
    KMfun=stats::stepfun(resKM$time,c(1,resKM$surv))
    resKMsurv=KMfun(Tobs_ord)
  } else {resKMsurv=resKM$surv}
  RMSTbase<-Rmst_ord(Tobs_ord,status_ord,tau)
  #Compute the RSMT on the vector of observed times (that is for tau=T_i, i=1, ..., n)
  integ=c(resKMsurv[1:(n-1)])*diff(Tobs_ord) #integ of S over [T_{i-1},T_i], i=2,...,n
  integcum=cumsum(integ) #integ of S over [T_1,T_i], i=2,...,n
  RMSTvect=c(Tobs_ord[1],Tobs_ord[1]+integcum) #RMST for all observations (censored and non censored)!!

  integ1=RMSTbase-RMSTvect #integ of S over Tl and tau
  index1=which(integ1>=0)#correspond to T_l<tau #(integ1>0)
  index2=which(integ1<0)#correspond to T_l>tau #rep(1,n)-index1#
  integ1[index2]<-0

  H<-seq(1,1/n,by=-1/n)
  part1=cumsum(status_ord/H^2)/n #integ of dH_1/(H^2) over [0,T_l]

  #construct integ of S over [T_i,T_l \wedge tau]
  n1=length(index1)
  n2=length(index2)
  shortintegcum=integcum[index1][1:(n1-1)] #integ of S over [T_1,T_i], i=2,...,n but for T-i<tau

  integ2tau=matrix(rep(shortintegcum,n1-1),ncol=n1-1)
  integ2tau=integ2tau-matrix(rep(c(0,shortintegcum[1:(n1-2)]),n1-1),byrow=TRUE,ncol=n1-1) #(n1-1)*(n1-1) matrix
  if (n2!=0){
  integ2tau=rbind(integ2tau,matrix(rep(integ1[1:(n1-1)],n2),byrow=TRUE,ncol=n1-1)) #(n-1)*(n1-1) matrix
  }
  integ3=integ2tau*matrix(rep(status_ord[1:(n1-1)]/(H[1:(n1-1)])^2,n-1)/n,byrow=TRUE,ncol=n1-1)
  integ3[upper.tri(integ3)] <- 0
  part2=apply(integ3,1,sum)
  if (n2!=0){
  part2[n1:(n-1)]<-part2[n1:(n-1)]+rep((integ1[n1]*status_ord[n1]/(H[n1])^2)/n,n-n1)
  }
  part2=c(0,part2)

  part3=integ1*status_ord/H

  VonM<-RMSTbase+part1*integ1+part2-part3
  return(VonM)

}


