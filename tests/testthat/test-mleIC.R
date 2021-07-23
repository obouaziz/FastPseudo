#Those tests are for the function mleIC but also for the pseudoIC function
#since this latter function needs an object rendered from mleIC as an input.

n=50
Left=runif(n,0,5)
Right=Left+runif(n,0,2)

test_that("mleIC throws an error", {
  expect_error(mleIC(Left,Right,cuts=c(1,1,2,3,4),a=rep(log(0.5),length(cuts)+1)))
  expect_error(mleIC(Left,c(Right,10),cuts=c(1,2,3,4),a=rep(log(0.5),5)))
  expect_error(mleIC(Left,Right,cuts=c(1,2,3,4),a=rep(log(0.5),6)))
})

##Only exact observations
#Exponential model
n=10000
Left<-Right<-rexp(n,0.5)
result=mleIC(Left,Right,cuts=NULL,a=rep(log(0.5),1),verbose=FALSE)
test_that("mleIC works with the exponential when there are only exact observations", {
  expect_equal(round(result$lambda,1),0.5)
})
#pch model
alpha=c(0.2,0.4,0.5,0.6)
Left<-Right<-rsurv(n,cuts=c(2,4,5),alpha=alpha)
result=mleIC(Left,Right,cuts=c(2,4,5),a=rep(log(0.5),4),verbose=FALSE)
test_that("mleIC works with the exponential when there are only exact observations", {
  expect_equal(round(max(abs(result$lambda-alpha)),1),0)
})

##Only interval-censored observations
n=4000
cuts=c(20,40,50)
alpha=c(0,0.025,0.05,0.1)
TrueTime=rsurv(n,cuts,alpha) #generate true data from the pch model
#Simulation of interval-censored data
Right<-rep(Inf,n)
nb.visit=5
visTime=0;visit=matrix(0,n,nb.visit+1)
visit=cbind(visit,rep(Inf,n))
visit[,2]=visit[,1]+stats::runif(n,0,20)#runif(n,0,5)
schedule=12
for (i in 3:(nb.visit+1))
{
  visit[,i]=visit[,i-1]+stats::runif(n,0,schedule*2)
}
Left<-visit[,(nb.visit+1)]
J=sapply(1:(n),function(i)cut(TrueTime[i],breaks=c(visit[1:(n),][i,]),
                              labels=1:(nb.visit+1),right=FALSE)) #sum(is.na(J)) check!
Left[1:(n)]=sapply(1:(n),function(i)visit[1:(n),][i,J[i]])
Right[1:(n)]=sapply(1:(n),function(i)visit[1:(n),][i,as.numeric(J[i])+1])
#View(data.frame(Left,Right,TrueTime)) #To see the generated data
#mean(Right==Inf) #percentage of right-censored data
#mean(Left<20 & Right!=Inf)
#mean(Left>50)
result=mleIC(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
             maxiter=1000,tol=1e-12,verbose=FALSE)
test_that("mleIC works when there are only interval censored observations", {
  expect_equal(1*(max(round(abs(result$lambda-alpha),2))<0.02),1)
})

test_that("pseudoIC should return a warning message",{
  expect_warning(pseudoIC(result,Left,Right,tseq=30))
})

#exponential model
#check if pseudoIC works for the exponential model
result=mleIC(Left,Right,cuts=NULL,a=rep(log(0.5),1),
             maxiter=1000,tol=1e-12,verbose=FALSE)
test_that("pseudoIC for the exponential model",{
  expect_error(pseudoIC(result,Left,Right,tseq=30),NA)
  pseudo_val=pseudoIC(result,Left,Right,tseq=30)
  expect_equal(mean(pseudo_val),exp(-pchcumhaz(30,cuts=NULL,result$lambda)))
  expect_error(pseudoIC(result,Left,Right,tseq=c(20,30,40)),NA)
})

test_that("pseudoIC with Left and Right not of same length for the exponential model",{
  expect_error(pseudoIC(result,Left,Right[-1],tseq=30))
})


#pch model
#check if pseudoIC works for the pch model
result=mleIC(Left,Right,cuts=c(40,50),a=rep(log(0.5),3),
             maxiter=1000,tol=1e-12,verbose=FALSE)
test_that("pseudoIC for the pch model",{
  expect_error(pseudoIC(result,Left,Right,tseq=30),NA)
  expect_error(pseudoIC(result,Left,Right,tseq=c(20,30,40)),NA)
  pseudo_val=pseudoIC(result,Left,Right,tseq=30)
  expect_equal(mean(pseudo_val),exp(-pchcumhaz(30,cuts=c(40,50),result$lambda)))
  expect_error(pseudoIC(result,Left,Right,tseq=c(20,30,40)),NA)
})

test_that("pseudoIC with Left and Right not of same length for the pch model",{
  expect_error(pseudoIC(result,Left,Right[-1],tseq=30))
})


##Mixed case of interval-censored and exact observations
Left[1:100]<-Right[1:100]<-TrueTime[1:100] #add 100 exact observations
result=mleIC(Left,Right,cuts=cuts,a=rep(log(0.5),length(cuts)+1),
             maxiter=1000,tol=1e-12,verbose=FALSE)
test_that("mleIC works when there are exact and interval censored observations", {
  expect_equal(1*(max(round(abs(result$lambda-alpha),2))<0.02),1)
})


