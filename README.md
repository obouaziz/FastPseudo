
# FastPseudo

<!-- badges: start -->
<!-- badges: end -->

Implementation of pseudo-values for right-censored or interval-censored using 
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
of cut values. It can also includes exact observations. The algorithm is based on the EM algorithm by considering the true 
times as unobserved variables.
- `pchcumhaz` computes the cumulative hazard function of a pch model. It takes as input a sequence of time at which the 
cumulative hazard needs to be computed, the cuts and the values of the hazard between two cuts.
- `pseudoIC` implements the pseudo-values for the survival function or the RMST. The fast approximation is based on the pch model. It takes as input interval-censored data and a value of the endpoint if the RMST should be computed.

## Installation

You can install the released version of FastPseudo with:

``` r
install_github("obouaziz/FastPseudo")
library(FastPseudo)
```

## Statistical Methodology

The programs are based on the following paper: [CRAN](https://CRAN.R-project.org)

## Examples

This is a basic example which shows you how to solve a common problem:

``` r
library(FastPseudo)
## basic example code
```

