# mix2mgsem: Double-Mixture Multigroup SEM

<!-- badges: start -->

[![R-CMD-check](https://github.com/hzhaoR/mix2mgsem/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hzhaoR/mix2mgsem/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`mix2mgsem` is an R package for two-step double-mixture multigroup structural equation modeling.

The method is designed for comparing structural relations across many groups while accounting for measurement model heterogeneity. It combines measurement model clustering and structural relations clustering in a two-step workflow. For a detailed description of the method, see: [MixMG-SEM with double mixture modeling](https://www.sciencedirect.com/science/article/pii/S2590260126000214).

## Overview

The current implementation consists of two steps:

1.  **Step 1: measurement model clustering**\
    Groups are clustered based on measurement model parameters using mixture multigroup factor analysis via `mixmgfa`.

2.  **Step 2: structural relations clustering**\
    Group-specific factor covariance matrices from Step 1 are used as input for clustering groups based on structural relations.

## Main functions

The package currently provides two main functions:

- `MixMix_Step1()`: runs measurement model clustering and extracts group-specific factor covariance matrices.
- `MixMix_Step2()`: runs structural relations clustering based on the Step 1 output.

The package also contains internal helper functions used for the EM algorithm and likelihood calculations.

## Schematic usage

``` r
library(mix2mgsem)

# Step 1: measurement model clustering
step1_model <- "
  F1 =~ x1 + x2 + x3 + x4 + x5
  F2 =~ z1 + z2 + z3 + z4 + z5
  F3 =~ m1 + m2 + m3 + m4 + m5
  F4 =~ y1 + y2 + y3 + y4 + y5
"

step1_out <- MixMix_Step1(
  data = dat,
  step1model = step1_model,
  group = "country",
  MM.cluster.spec = c("loadings"),
  MM.nclus = 2,
  MM.design = design,
  invar_loadings = invar_loadings,
  markers = markers
)

# Step 2: structural relations clustering
step2_model <- "
  F4 ~ F1 + F3
  F3 ~ F1 + F2
"

step2_out <- MixMix_Step2(
  s1out = step1_out,
  step2model = step2_model,
  nclus = 2,
  nfree = 4,
  seed = 100
)
```

## Development status

`mix2mgsem` is currently under active development. The current version includes:

- an R package structure with `DESCRIPTION`, `NAMESPACE`, and `man/` documentation files;
- roxygen2 documentation for the main workflow functions;
- basic `testthat` unit tests for internal E-step behavior and input validation;
- GitHub Actions continuous integration using R CMD check.

Planned improvements include:

- a small runnable example;
- additional tests for expected output structure;
- smoke tests for `MixMix_Step2()` using a minimal mock Step 1 object;
- simulation-based validation examples;
- further refactoring of the Step 1 and Step 2 code.

## Installation

The package is not yet on CRAN. The development version can be installed from GitHub with:

``` r
# install.packages("devtools") 
devtools::install_github("hzhaoR/mix2mgsem") 
```

## Acknowledgements

`mix2mgsem` builds upon several existing R packages, including:

- **mixmgfa** for mixture multigroup factor analysis;
- **mmgsem** for mixture multigroup SEM;
- **lavaan** for structural equation modeling.

Parts of the Step 2 implementation are adapted from the `mmgsem` package developed by Andres Felipe Perez Alonso and extended for the double-mixture multigroup SEM implemented in `mix2mgsem`.

The Step 1 scaling/rotation helper is adapted from updated \`mixmgfa\` code provided by Kim De Roover and included internally to support the current two-step workflow.

## Funding

This package was developed as part of the ERC-funded project (PROCESSHETEROGENEITY, 101040754, awarded to Kim De Roover). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.
