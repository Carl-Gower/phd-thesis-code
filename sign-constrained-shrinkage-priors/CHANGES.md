# Change summary

## Common interface

The modified samplers accept `pos_con` and `neg_con` vectors of one-based coefficient indices. The release preparation added checks for whole-number, in-range, non-overlapping indices and preserved the original unconstrained route when both are empty.

## R2D2 and Dirichlet-Laplace

`dl()` uses coordinate-wise truncated-normal coefficient updates when signs are supplied. `r2d2marg()` draws from a truncated multivariate normal. The starting coefficient vector is placed inside the requested support so that all returned rows obey the constraints. The overlay also records the added package dependencies.

## bayesreg

The scale-mixture coefficient update uses a truncated multivariate normal for constrained fits. Unsupported Metropolis-Hastings response models stop with an explanatory error. A fallback to one core was added for systems on which R cannot detect the available processor count.

## Boom and BoomSpikeSlab

The R interface validates and passes the constraints to `BregVsSampler`. The classic Gaussian SSVS coefficient update uses a coordinate-wise truncated-normal sweep over the included coefficients. Adaptive SSVS, ODA and Student-error calls with constraints are rejected because those sampler classes have no corresponding implementation.

## Normal-Gamma

`NG.m` retains the original block draw when unconstrained and uses a coordinate-wise truncated-normal update when signs are supplied. Input validation was added, the truncated-normal helper now uses MATLAB's distribution truncation, and the function no longer changes the global warning state.

## Installation

The installer downloads verified upstream sources, applies the overlays and installs the modified packages in a separate local R library. Modified versions carry a `.9000` suffix. `use_overlays.R` selects that library and checks for packages already loaded from another location. The installation checks call the installed package functions, including the compiled spike-and-slab sampler.

The Normal-Gamma implementation includes the unmodified RANDRAW dependency and its BSD licence in `normal-gamma/private/`.
