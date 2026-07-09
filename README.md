# mix2mgsem: Double-Mixture Multigroup SEM

<!-- badges: start -->

[![R-CMD-check](https://github.com/hzhaoR/mix2mgsem/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hzhaoR/mix2mgsem/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`mix2mgsem` is an R package for two-step double-mixture multigroup structural equation modeling.

The method is designed for comparing structural relations across many groups while accounting for measurement non-invariance. It combines measurement model clustering and structural relations clustering in a two-step workflow. For a detailed description of the method, see: [MixMG-SEM with double mixture modeling](https://www.sciencedirect.com/science/article/pii/S2590260126000214).

## Overview

The current implementation consists of two steps:

1.  **Step 1: measurement model clustering**\
    Groups are clustered based on (a subset of) measurement model (MM) parameters using mixture multigroup factor analysis via `mixmgfa`.

2.  **Step 2: structural relations clustering**\
    Group-specific factor covariance matrices from Step 1 are used as input for clustering groups based on structural relations (SR).

## Installation

The package is not yet on CRAN. The development version can be installed from GitHub with:

``` r
# install.packages("devtools") 
devtools::install_github("hzhaoR/mix2mgsem") 
```

## Example

This example uses simulated data included with the package. It runs the two-step workflow: first MM-clustering with `MixMix_Step1()`, then SR-clustering with `MixMix_Step2()`.

``` r
library(mix2mgsem)

# load data
sim <- readRDS(system.file(
  "extdata",
  "data19346.rds",
  package = "mix2mgsem"
))

dat <- as.data.frame(sim$SimData)

# step 1 model
step1model <- '
  F1 =~ x1 + x2 + x3 + x4 + x5
  F2 =~ z1 + z2 + z3 + z4 + z5
  F3 =~ m1 + m2 + m3 + m4 + m5
  F4 =~ y1 + y2 + y3 + y4 + y5
'

# loading matrix and marker variables 
loading_matrix <- matrix(
  c(rep(c(rep(1, 5), rep(0, 20)), 3), rep(1, 5)),
  nrow = 20,
  ncol = 4,
  dimnames = list(
    c(paste0("x", 1:5),
      paste0("z", 1:5),
      paste0("m", 1:5),
      paste0("y", 1:5)),
    paste0("F", 1:4)
  )
)
markers <- loading_matrix
markers[2:5, 1] <- markers[7:10, 2] <- markers[12:15, 3] <- markers[17:20, 4] <- 0

# run step 1 with 2 MM-clusters
s1out <- MixMix_Step1(
  data = dat,
  step1model = step1model,
  group = "group",
  MM.cluster.spec = "loadings",
  MM.nclus = 2,
  MM.design = loading_matrix,
  markers = markers,
  seed = 100
)

# step 2 model
step2model <- '
  F4 ~ F1 + F3
  F3 ~ F1 + F2
'

# run step 2 with 2 SR-clusters
s2out <- MixMix_Step2(
  s1out = s1out,
  step2model = step2model,
  nclus = 2,
  nfree = 16,
  seed = 100,
  printing = FALSE
)

# s2out$posteriors
```

## Main functions

The package currently provides two main functions:

- `MixMix_Step1()`: runs MM-clustering and extracts group-specific factor covariance matrices.
- `MixMix_Step2()`: runs SR-clustering based on the Step 1 output.

The package also contains internal helper functions used for the EM algorithm and likelihood calculations.

## Development status

`mix2mgsem` is currently under active development. The current version includes:

- validation tests for both steps using simulated data, checking recovery of MM-clusters, factor loadings, SR-clusters, and cluster-specific regression coefficients;

Planned improvements include:

- further refactoring of the Step 1 and Step 2 code.

## Acknowledgements

`mix2mgsem` builds upon several existing R packages, including:

- **mixmgfa** for mixture multigroup factor analysis;
- **mmgsem** for mixture multigroup SEM;
- **lavaan** for structural equation modeling.

Parts of the Step 2 implementation are adapted from the `mmgsem` package developed by Andres Felipe Perez Alonso and extended for the double-mixture multigroup SEM implemented in `mix2mgsem`.

The Step 1 scaling/rotation helper is adapted from updated `mixmgfa` code provided by Kim De Roover and included internally to support the current two-step workflow. `MixMix_Step1()` can be run with the public version of `mixmgfa` when `invar_loadings = NULL`; or it requires a development version of `mixmgfa` that supports that argument.

## Funding

This package was developed as part of the ERC-funded project (PROCESSHETEROGENEITY, 101040754, awarded to Kim De Roover). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.
