library(survival)
n=1000
lambda=0.5
cpar=0.2
tau=2
#set.seed(30)
TrueTime=rexp(n,lambda)
Cens=rexp(n,cpar)

Tobs=pmin(TrueTime,Cens)
Tsort<-sort(Tobs,index.return=TRUE)
Tobs_ord<-Tsort$x

status=TrueTime<=Cens #mean(status) #0.7
status_ord<-status[Tsort$ix]

our_estimate=Rmst(Tobs,status,tau)
our_estimate2=Rmst_ord(Tobs_ord,status_ord,tau)
km=survfit(Surv(Tobs,status)~1)
#print(survfit(Surv(Tobs,status)~1),print.rmean=TRUE, rmean=2)
survival_rmst=survival:::survmean(km, rmean=tau)[[1]]["*rmean"]
truth=(1-exp(-lambda*tau))/lambda

test_that("Rmst function works", {
  expect_equal(our_estimate$tau,tau)
  expect_equal(our_estimate$rmst,survival_rmst[[1]])#round(our_estimate$rmst,2),round(survival_rmst[[1]],2)
  expect_equal(round(our_estimate2,20),round(survival_rmst[[1]],20))
  expect_equal(sum(abs(our_estimate$rmst-truth)<2),1)
})





