
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

Here is the output (we only display the effect of X). We obtain an effect very 
close to the true value.

Call:
geese(formula = Y ~ X + as.factor(Time) - 1, id = id, data = data_pseudo1, 
    family = "gaussian", mean.link = "cloglog", scale.fix = TRUE, 
    corstr = "independence", jack = TRUE)

Mean Model:
 Mean Link:                 cloglog 
 Variance to Mean Relation: gaussian 

 Coefficients:

|               | estimate    |     san.se |     ajs.se  |   wald   |        p
| ------------- |:----:|:----:|:----:|:----:|:----:| 
X               | -0.70061472 | 0.07440208 | 0.07434210  | 88.672378| 0.000000e+00 


Scale is fixed.

Correlation Model:
 Correlation Structure:     independence 

Returned Error Value:    0 
Number of clusters:   4000   Maximum cluster size: 10 

We now illustrate the pseudo-values for the RMST. We first simulate data according to 
a linear model with two covariates. It can be shown that the corresponding RMST will 
then follow a linear relationship with respect to the two covariates and their 
two interactions. The true effects were empirically estimated on a sample of size 
1e7 and were found to be equal to 3.812552, 0.05705492, 0.05730445, 0.1044318.

``` r
#Simulation in a linear model
set.seed(28)
n<-10000
sigma=3
cpar=0.07
tau=4
alpha=c(5.5,0.25,0.25)
X1=rbinom(n,1,0.5)
X2=rbinom(n,1,0.5)
epsi=runif(n,-sigma,sigma)
TrueTime=alpha[1]+alpha[2]*X1+alpha[3]*X2+epsi
Cens=rexp(n,cpar)
Tobs=pmin(TrueTime,Cens)
Tsort<-sort(Tobs,index.return=TRUE)
Tobs_ord<-Tsort$x
status=TrueTime<=Cens #approximately 35% of censoring
status_ord<-status[Tsort$ix]
X_ord1=X1[Tsort$ix];X_ord2=X2[Tsort$ix]
X00=(X1==0 & X2==0)
X01=(X1==0 & X2==1)
X10=(X1==1 & X2==0)
X11=(X1==1 & X2==1)
X_ord00=X00[Tsort$ix];X_ord01=X01[Tsort$ix];X_ord10=X10[Tsort$ix];X_ord11=X11[Tsort$ix]
#The true value of the parameters are c(3.812552,0.05705492,0.05730445,0.1044318)

pseudo_val=pseudoKM(Tobs,status,tau)$pseudoval
data_pseudo<-data.frame(Y=c(pseudo_val),X2=X01,X3=X10,X4=X11,id=rep(1:n))
result=geese(Y~X2+X3+X4,id=id,jack=TRUE,family="gaussian",
mean.link="identity", corstr="independence",scale.fix=TRUE,data=data_pseudo)
summary(result)
```

Call:
geese(formula = Y ~ X2 + X3 + X4, id = id, data = data_pseudo, 
    family = "gaussian", mean.link = "identity", scale.fix = TRUE, 
    corstr = "independence", jack = TRUE)

Mean Model:
 Mean Link:                 identity 
 Variance to Mean Relation: gaussian 

 Coefficients:
 
 |               | estimate    |     san.se |     ajs.se  |   wald   |        p
| ------------- |:----:|:----:|:----:|:----:|:----:| 
(Intercept) | 3.82269101 | 0.008294085 | 0.008295309 | 212422.96376 | 0.0000000000
X2TRUE      | 0.04045047 | 0.010693009 | 0.010694562     | 14.31025 | 0.0001550181
X3TRUE      | 0.04017529 | 0.010816528 | 0.010818171     | 13.79565 | 0.0002038073
X4TRUE      | 0.09621505 | 0.009633681 | 0.009635118     | 99.74738 | 0.0000000000

Scale is fixed.

Correlation Model:
 Correlation Structure:     independence 

Returned Error Value:    0 
Number of clusters:   10000   Maximum cluster size: 1 



<!-- badges: end -->
