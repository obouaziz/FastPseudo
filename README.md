
# FastPseudo

<!-- badges: start -->

Implementation of pseudo-values for right-censored or interval-censored data using 
a fast (and very accurate) approximation. Both the pseudo-values for the survival function
and for the Restricted Mean Survival Time (RMST) have been implemented. For right-censored 
observations the algorithms are based on the Kaplan-Meier estimator while for interval-censored 
data they are based on the piecewise constant hazard (PCH) model.

The following functions are for right-censored observations only: 

- `pseudoKM` implements the pseudo-values for the Kaplan-Meier estimator. It can returns the pseudo-values for 
the survival function or for the RMST. It takes as input a continuous time variable, its censoring indicator 
and a value of the endpoint if the RMST should be computed.
- `Rmst` computes the RMST for the Kaplan-Meier estimator. It takes as input a continuous time variable, its censoring indicator
and a value of the endpoint for the RMST. 

The following functions are for interval-censored observations (which can contain right-censored data):
- `rsurv` simulates data following a pch model. Only the cuts and values of the hazard between two cuts must be specified.
- `RmstIC` computes the RMST based on the pch model. Only the cuts, values of the hazard between two cuts and 
a value of the endpoint for the RMST must be specified.
- `mleIC` computes the maximum likelihood estimator from a pch model. It takes as input interval-censored data and a sequence 
of cut values. It can also include exact observations. The estimation method is based on the EM algorithm by considering the true 
times as unobserved variables.
- `pchcumhaz` computes the cumulative hazard function of a pch model. It takes as input a sequence of times at which the 
cumulative hazard needs to be computed, the cuts and the values of the hazard between two cuts.
- `pseudoIC` implements the pseudo-values for the survival function or the RMST. The fast approximation is based on the pch model. It takes as input interval-censored data and a value of the endpoint if the RMST should be computed.

## Installation

You can install the released version of FastPseudo, from the `devtools` package with:

``` r
install_github("obouaziz/FastPseudo")
library(FastPseudo)
```

## Statistical Methodology

The programs are based on the following paper: [CRAN](https://CRAN.R-project.org). 
While the paper mention that the pseudo-values are approximations of the standard 
jackknife method, it is important to stress that the approximation is extremely 
accurate even for small sample sizes.

## Examples for right-censored data

We start with a very simple example. We first true survival times following the Weibull 
distribution and a censoring variable giving approximately $21\%$ of censoring. We 
then use the function `pseudoKM` to compute pseudo values for the survival function. Plots for 
for two non-censored pseudo values and two censored pseudo values are displayed.

``` r
require(devtools)
install_github("obouaziz/FastPseudo")
library(FastPseudo)
require(survival)
n=100
cpar=0.1
set.seed(28)
TrueTime=rweibull(n,shape=0.5,scale=2)
Cens=rexp(n,cpar)
Tobs=pmin(TrueTime,Cens) #observed times
status=TrueTime<=Cens #mean(status) #21% of censoring on average

#use tau=NULL to specify that the pseudo-values are for the survival
#function and not for the RMST.
pseudo_val=pseudoKM(Tobs,status,tau=NULL)
pseudoval=pseudo_val$pseudoval
tseq=pseudo_val$tseq

par(mfrow=c(2,2))
plot(tseq,pseudoval[,3],type="s",xlab="Time",ylab="Pseudo-val. for obs. 3")
abline(v=Tobs[3],lty=2,col="red")
plot(tseq,pseudoval[,48],type="s",xlab="Time",ylab="Pseudo-val. for obs. 48")
abline(v=Tobs[48],lty=2,col="red")
plot(tseq,pseudoval[,1],type="s",xlab="Time",ylab="Pseudo-val. for obs. 1")
abline(v=Tobs[1],lty=2,col="red")
plot(tseq,pseudoval[,5],type="s",xlab="Time",ylab="Pseudo-val. for obs. 5")
abline(v=Tobs[5],lty=2,col="red")
par(mfrow=c(1,1))
```

![](Image/pseudoObs.png)


We now simulate right-censored data from a Cox model with Weibull baseline. There 
is only one covariate whose effect is equal to -0.69315. We compute the pseudo-values 
at 10 different time points and we implement generalised estimating equations using the 
`geepack`package with link function equal to log(-log(.)). We should retrieve the 
value of the parameter.
``` r
set.seed(28)
require(geepack)
require(survival)
n=4000
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
SurvEst<-matrix(NA,n,n)

#Compute pseudo-values using pseudoKM
pseudo=pseudoKM(Tobs,status,tau=NULL)
#the pseudo values. Individuals are on the columns, while the time is on the rows.
pseudo_val=pseudo$pseudoval 
pseudo_val=pseudo_val[n*seq(5,95,by=10)/100,] #we compute the pseudo-values for 10 times
tseq=pseudo$tseq[n*seq(5,95,by=10)/100]
M=10
data_pseudo1<-data.frame(Y=1-c(pseudo_val),X=rep(X,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
#data_pseudo1<-data.frame(Y=1-c(pseudo_val),X=rep(X_ord,each=M),Time=rep(tseq,n),id=rep(1:n,each=M))
result<-geese(Y~X+as.factor(Time)-1,id=id,jack=TRUE,family="gaussian",
mean.link="cloglog", corstr="independence",scale.fix=TRUE,data=data_pseudo1)
summary(result)
```

Call:
geese(formula = Y ~ X + as.factor(Time) - 1, id = id, data = data_pseudo1, 
    family = "gaussian", mean.link = "cloglog", scale.fix = TRUE, 
    corstr = "independence", jack = TRUE)

Mean Model:
 Mean Link:                 cloglog 
 Variance to Mean Relation: gaussian 

 Coefficients:
                                   estimate     san.se     ajs.se        wald            p
X                               -0.70061472 0.07440208 0.07434210   88.672378 0.000000e+00
as.factor(Time)4.58335382087617 -6.53429441 0.50775420 0.50715275  165.611446 0.000000e+00
as.factor(Time)11.8286543020646 -2.74044302 0.08585274 0.08575709 1018.903980 0.000000e+00
as.factor(Time)16.6755085579274 -1.78659252 0.06120366 0.06114036  852.111025 0.000000e+00
as.factor(Time)20.5772930286843 -1.16647825 0.05141861 0.05136871  514.650789 0.000000e+00
as.factor(Time)24.0869086570905 -0.70279219 0.04624723 0.04620462  230.931094 0.000000e+00
as.factor(Time)27.3635276924689 -0.30775626 0.04378098 0.04374176   49.413186 2.073453e-12
as.factor(Time)30.8475008520177  0.04348487 0.04243428 0.04239671    1.050129 3.054774e-01
as.factor(Time)34.3261917181898  0.39800295 0.04198368 0.04194640   89.869366 0.000000e+00
as.factor(Time)38.7245546923862  0.77791031 0.04340787 0.04336864  321.160505 0.000000e+00
as.factor(Time)46.0320710204542  1.29594973 0.04960845 0.04956259  682.440899 0.000000e+00

Scale is fixed.

Correlation Model:
 Correlation Structure:     independence 

Returned Error Value:    0 
Number of clusters:   4000   Maximum cluster size: 10 


<!-- badges: end -->
