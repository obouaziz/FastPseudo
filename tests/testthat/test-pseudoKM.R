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

VonM=pseudoKM(Tobs_ord,status_ord,tau=NULL)
pseudo_VM=VonM[n*seq(5,95,by=10)/100,] #we compute the pseudo-values for 10 times
tseq=Tobs_ord[n*seq(5,95,by=10)/100]
M=length(n*seq(5,95,by=10)/100)
data_pseudo<-data.frame(Y=1-c(pseudo_VM),X=rep(X_ord,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
result<-geese(Y~X+as.factor(Time)-1,id=id,jack=TRUE,family="gaussian", mean.link="cloglog", corstr="independence",scale.fix=TRUE,data=data_pseudo)

test_that("pseudoKM for survival function works", {
  expect_equal(round(result$beta[1][[1]],2),-0.55)
})
