# mix2mgsem: Double-Mixture Multigroup SEM

`mix2mgsem` is an R package for two-step double-mixture multigroup structural equation modeling.

The method is designed for comparing structural relations across many groups while accounting for measurement-model heterogeneity. It combines measurement-model clustering and structural-relation clustering in a two-step workflow. For a detailed description of the method, see: [MixMG-SEM with double mixture modeling](https://www.sciencedirect.com/science/article/pii/S2590260126000214).

## Overview

The current implementation consists of two steps:

1.  **Step 1: measurement-model clustering**\
    Groups are clustered based on measurement model parameters using mixture multigroup factor analysis via `mixmgfa`.

2.  **Step 2: structural-relation clustering**\
    Group-specific factor covariance matrices from Step 1 are used as input for clustering groups based on structural relations.

This workflow is useful when researchers want to compare latent structural relations across many groups, but measurement parameters may differ across groups or measurement clusters.

## Development status

This package is currently under development.

At the moment, the package contains early implementations of:

- `MixMix_Step1()`: measurement model clustering and group-specific factor covariance matrices.
- `MixMix_Step2()`: structural relations clustering based on the Step 1 output.

The interface, documentation, tests, and examples are still being developed.

## Installation

The package is not yet on CRAN. Once the GitHub version is ready, it can be installed with:

``` r
# install.packages("devtools") 
devtools::install_github("hzhaoR/mix2mgsem") 
```

## Acknowledgements

`mix2mgsem` builds upon several existing R packages, including:

- **mixmgfa** for mixture multigroup factor analysis.
- **mmgsem** for mixture multigroup SEM.
- **lavaan** for structural equation modeling.

Parts of the Step 2 implementation are adapted from the `mmgsem` package developed by Andres Felipe Perez Alonso and extended for the double-mixture multigroup SEM implemented in `mix2mgsem`.

## Funding

This package was developed as part of the ERC-funded project (PROCESSHETEROGENEITY, 101040754, awarded to Kim De Roover). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.
