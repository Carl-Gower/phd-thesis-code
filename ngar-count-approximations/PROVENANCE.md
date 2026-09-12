# Provenance

## Original method and code

The sampler implements the normal-gamma autoregressive dynamic regression model from Kalli and Griffin (2014). Carl Gower obtained the accompanying `NGAR.m` and `kf_NGAR.m` MATLAB files from Jim Griffin's former University of Kent webpage. The source copy inspected while preparing this repository did not contain a software licence.

The original files are not included separately in this folder. The uploaded runtime files are derivatives of them and are published with the permission described below.

## Changes made by Carl Gower

The uploaded `src/NGAR.m` adds optional normal and Poisson approximations to latent count draws. These approximations provide an alternative when the ordinary draws become computationally impractical at large counts or Poisson means. Below the selected threshold, or when no approximation condition applies, the ordinary MATLAB route is retained.

The uploaded version also provides one interface for the ordinary and approximate routes, optional output storage, run diagnostics, input checks and numerical safeguards. [`CHANGES.md`](CHANGES.md) gives a direct comparison with the original implementation.

## Permission

Maria Kalli and Jim E. Griffin gave Carl Gower written permission by email to put the modified code on GitHub. The correspondence is retained privately. It did not assign a named standard software licence; [`LICENSE.md`](LICENSE.md) explains what this means for reuse.

## Reference

Kalli, M. and Griffin, J. E. (2014). Time-varying sparsity in dynamic regression models. *Journal of Econometrics*, 178, 779–793. <https://doi.org/10.1016/j.jeconom.2013.10.012>
