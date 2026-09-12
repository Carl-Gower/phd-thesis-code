# PhD thesis code

Research code accompanying my PhD thesis, *Bayesian Shrinkage Priors for Time Series Analysis* (2025). It contains three sets of Bayesian regression methods that I adapted or extended during the thesis.

## Contents

| Folder | Language | What is included |
| --- | --- | --- |
| [`dsp-extensions/`](dsp-extensions/) | R | A standalone implementation based on Kowal, Matteson and Ruppert's dynamic shrinkage process code, with additional state-evolution and persistence options. |
| [`ngar-count-approximations/`](ngar-count-approximations/) | MATLAB | An implementation based on Kalli and Griffin's NGAR sampler, including the optional count approximations developed for my thesis. |
| [`sign-constrained-shrinkage-priors/`](sign-constrained-shrinkage-priors/) | R, C++ and MATLAB | Modifications that allow selected regression coefficients to be restricted to positive or negative values under several static shrinkage priors. |

To run a method, open its folder above and follow the README. Each README gives the requirements, a synthetic example, tests and an account of the changes to the original implementation.

## Scope

The examples demonstrate the reusable methods on synthetic data. Thesis datasets, saved results and the scripts for individual tables and figures are outside the scope of this release.

## Attribution and licences

Substantial parts of the implementations began with code written by the authors cited in the component READMEs. My contributions are the extensions and implementation changes described in each component’s change summary.

The repository has no single blanket licence because the three sections have different origins and licence terms. See [`LICENSE.md`](LICENSE.md) and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) before reusing or redistributing anything.

The NGAR MATLAB sources are published here with written permission from Maria Kalli and Jim E. Griffin. The Normal-Gamma MATLAB source is published with written permission from Jim E. Griffin. No named standard software licence was assigned to those components; their permission notices record the basis on which they appear here.

## Citation

If you use one of the methods, please cite both the relevant original paper and this code release. The paper references are given in each component folder; repository citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Contact

Questions and reproducibility reports can be opened as a GitHub issue.
