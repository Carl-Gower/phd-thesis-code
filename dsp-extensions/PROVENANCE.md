# Source provenance

## Upstream implementation

This component is derived from version 0.1.0 of Daniel R. Kowal's [`dsp`](https://github.com/drkowal/dsp) R package, distributed under GPL-2. The implementation accompanies:

> Kowal, D. R., Matteson, D. S. and Ruppert, D. (2019). Dynamic shrinkage processes. *Journal of the Royal Statistical Society: Series B*, 81(4), 781–804. <https://doi.org/10.1111/rssb.12325>

The upstream reference point is commit [`d5ea2548cb251c078b79f57a5a9e273d2650f4bd`](https://github.com/drkowal/dsp/tree/d5ea2548cb251c078b79f57a5a9e273d2650f4bd). The relevant upstream files were compared with that commit while preparing this release:

- `R/component_samplers.R`
- `R/helper_functions.R`
- `R/mcmc_samplers.R`
- `DESCRIPTION`

## Relationship to the upstream code

`R/DSP_rho.R` is a standalone derivative rather than an official release of the package. It brings the functions needed for the univariate and time-varying regression samplers into one file and adds the AR/DAR and `rho` extensions described in [`CHANGES.md`](CHANGES.md).

Retained functions were compared with their upstream counterparts. The source also preserves two attributions present in `dsp`: `uni.slice()` is based on code supplied by Radford M. Neal, and `rig()` is based on code from the `mgcv` package. Those routines remain part of this GPL-2 derivative.

## Authorship

- Original `dsp` implementation: Daniel R. Kowal.
- Dynamic shrinkage process paper: Daniel R. Kowal, David S. Matteson and David Ruppert.
- AR/DAR specifications, shared or predictor-specific `rho`, and the standalone integration in this component: Carl Gower.

## Licence

The upstream package declares `License: GPL-2`. This modified component is consequently distributed under GPL-2.0-only; the full terms are in [`LICENSE`](LICENSE).
