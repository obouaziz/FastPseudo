set.seed(28)
require(geepack)
n=100
shape=3;scale=30;cpar=0.01
#Simulate covariates
X=runif(n,0,1)
#simulate Cox model with Weibull baseline distribution
TrueTime=scale*(rexp(n,1)/exp(log(0.5)*X))^(1/shape)
#True regression parameter is log(0.5)=-0.6931472
Cens=rexp(n,cpar)
Tobs=pmin(TrueTime,Cens) #observed times
Tsort<-sort(Tobs,index.return=TRUE)
Tobs_ord<-Tsort$x
status=TrueTime<=Cens #mean(status) #25% of censoring on average
status_ord<-status[Tsort$ix]
X_ord=X[Tsort$ix]

VonM=pseudoKM(Tobs_ord,status_ord,tau=NULL)$pseudoval
pseudo_VM=VonM[n*seq(5,95,by=10)/100,] #we compute the pseudo-values for 10 times
tseq=Tobs_ord[n*seq(5,95,by=10)/100]
M=length(n*seq(5,95,by=10)/100)
data_pseudo<-data.frame(Y=1-c(pseudo_VM),X=rep(X_ord,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
result<-geese(Y~X+as.factor(Time)-1,id=id,jack=TRUE,family="gaussian", mean.link="cloglog", corstr="independence",scale.fix=TRUE,data=data_pseudo)

test_that("pseudoKM for survival function works: check beta value on a regression model", {
  expect_equal(round(result$beta[1][[1]],2),-0.55)
})


n=100
cpar=0.1
TrueTime=rweibull(n,shape=0.5,scale=2)
Cens=rexp(n,cpar)
Tobs=pmin(TrueTime,Cens) #observed times
status=TrueTime<=Cens #mean(status) #21% of censoring on average

pseudo_val=pseudoKM(Tobs,status,tau=NULL)
VonM=pseudo_val$pseudoval
tseq=pseudo_val$tseq

test_that("pseudoKM for survival function works: check survival estimate", {
  expect_equal(round(max(abs(apply(VonM,1,mean)-survfit(Surv(Tobs,status)~1)$surv)),4),0)
})

n<-100
lambda=0.5
TrueTime=rexp(n,rate=lambda)
tau=10
tau=c(2,4,6,8,10)
cpar=0.2
Cens=rexp(n,cpar)
Tobs=pmin(TrueTime,Cens) #observed times
status=TrueTime<=Cens #mean(status) #36% of censoring on average

pseudo_val=sapply(tau,function(x){pseudoKM(Tobs,status,tau=x)$pseudoval})
rmst_pseudo=apply(pseudo_val,2,mean) #check that mean of pseudo-values return the RMST
Rmst_est=sapply(tau,function(x){Rmst(Tobs,status,x)$rmst}) #compute the RMST

test_that("pseudoKM for RMST works: check on RMST", {
  expect_equal(round(max(abs(rmst_pseudo-Rmst_est)),4),0)
})


require(geepack)
n<-400
sigma=1
tau=400
alpha=c(3)
X=rnorm(n,2,0.25)
epsi=rnorm(n,0,sigma)
TrueTime=alpha*X+epsi
Cens=rep(Inf,n)
Tobs=pmin(TrueTime,Cens)
status=TrueTime<=Cens

VonM=pseudoKM(Tobs,status,tau)$pseudoval
data_VonM<-data.frame(Y=c(VonM),X=X,id=rep(1:n))
resultEst=geese(Y~X,id=id,jack=TRUE,family="gaussian", mean.link="identity",
                corstr="independence",scale.fix=TRUE,data=data_VonM)
summary(resultEst)

resultOrac=lm(TrueTime~X)
summary(resultOrac)

test_that("pseudoKM for RMST works: check the regression estimate in a linear model", {
  expect_equal(1*(round(abs(resultEst$beta[2][[1]]-resultOrac$coefficients[2][[1]]),4)<0.2),1)
})

test_that("Time and status must have the same length",{
  expect_error(pseudoKM(c(1,2,3,4),c(1,1,0)))
})

test_that("Time and status must have the same length (tau is specified)",{
  expect_error(pseudoKM(c(1,2,3,4),c(1,1,0),tau=2.5))
})

test_that("Time and status must have the same length (tau is specified)",{
  expect_error(pseudoKM(c(1,2,3,4),c(1,1,0,0),tau=Inf))
})

test_that("Time and status must have the same length (tau is specified)",{
  expect_error(pseudoKM(c(1,2,3,4),c(1,1,0,1),tau=NA))
})

test_that("Time and status must have the same length (tau is specified)",{
  expect_error(pseudoKM(c(1,2,3,4),c(1,1,0,1),tau=2.5))
})

