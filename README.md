
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

[1] 1.264241 1.729329 1.900426 1.963369 1.986524

<!-- badges: end -->
