# Source provenance

The R overlays contain the modified source files and the metadata needed to install them. The MATLAB implementation also includes the licensed RANDRAW dependency. Compiled objects, thesis data and saved results are omitted.

| Component | Reference source | Files changed in this repository | Terms |
| --- | --- | --- | --- |
| `bayesreg` | CRAN version 1.3 | `DESCRIPTION`, `R/bayesreg.R` | GPL-3 or later |
| R2D2 | [`yandorazhang/R2D2`](https://github.com/yandorazhang/R2D2), commit `e734639929abb60e616c114ac7fe4e2beb5c7f9d` | `DESCRIPTION`, `R/dl.R`, `R/r2d2marg.R` | MIT |
| `Boom` | CRAN version 0.9.16 | `BregVsSampler.hpp`, `BregVsSampler.cpp` | LGPL-2.1 |
| `BoomSpikeSlab` | CRAN version 1.2.7 | `R/lm.spike.R`, `src/spike_slab_wrapper.cc` | LGPL-2.1 |
| Normal-Gamma | `NG.m` supplied privately by Jim E. Griffin | `NG.m` | Written permission to publish on GitHub; no named standard software licence assigned |

## Nature of the changes

- `bayesreg` and R2D2 add optional sign-index arguments and use truncated-normal coefficient updates when constraints are present. Their original unconstrained update paths are retained.
- The `Boom`/`BoomSpikeSlab` overlay passes sign indices from R to the classic SSVS sampler and applies coordinate-wise truncated-normal coefficient updates. Unsupported sampler routes stop rather than silently ignoring requested restrictions.
- `NG.m` retains the original unconstrained update and uses a coordinate-wise truncated-normal update when constraints are supplied. It also adds input checks and avoids changing MATLAB's global warning state.

The modifications were made by Carl Gower. They do not change the authorship of the upstream packages or methods.

## Normal-Gamma permission record

Jim E. Griffin gave Carl Gower written permission by email to put the modified code on GitHub. The correspondence is retained privately. It did not assign a named standard software licence, so the repository records the permission without interpreting broader terms. See [`normal-gamma/LICENSE.md`](normal-gamma/LICENSE.md).
