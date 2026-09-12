# Third-party code and provenance

This repository distinguishes the original implementations from the changes made by Carl Gower for his PhD. More detailed, file-level notes and change summaries are included in each component folder.

| Component | Original work | Local contribution | Source/version | Licence or permission |
| --- | --- | --- | --- | --- |
| Dynamic shrinkage processes | Daniel R. Kowal's `dsp` implementation accompanying Kowal, Matteson and Ruppert (2019) | Standalone integration, additional AR/DAR evolution choices, shared or predictor-specific persistence, and implementation safeguards | [`drkowal/dsp`](https://github.com/drkowal/dsp), commit `d5ea2548cb251c078b79f57a5a9e273d2650f4bd` | GPL-2 |
| NGAR | NGAR sampler by Maria Kalli and Jim E. Griffin accompanying Kalli and Griffin (2014) | Optional normal and Poisson count approximations, a unified interface, diagnostics and numerical safeguards | Copy originally obtained from Jim E. Griffin's University of Kent webpage | Written permission from Maria Kalli and Jim E. Griffin to publish on GitHub; no named standard software licence assigned |
| `bayesreg` | Enes Makalic and Daniel F. Schmidt | Truncated-normal coefficient updates for requested signs | `bayesreg` 1.3 | GPL (>= 3) |
| `Boom` | Steven L. Scott and contributors | Sign restrictions in the relevant spike-and-slab sampler | `Boom` 0.9.16 | LGPL-2.1 |
| `BoomSpikeSlab` | Steven L. Scott and contributors | R/C++ interface changes for sign restrictions | `BoomSpikeSlab` 1.2.7 | LGPL-2.1 |
| R2-D2 code | Yan Dora Zhang and co-authors | Sign-constrained variants of selected samplers | [`yandorazhang/R2D2`](https://github.com/yandorazhang/R2D2), commit `e734639929abb60e616c114ac7fe4e2beb5c7f9d` | MIT |
| Normal-Gamma | Implementation supplied privately by Jim E. Griffin, associated with Griffin and Brown (2010) | Truncated-normal coefficient update for requested signs | Emailed source copy; no public version identifier | Written permission from Jim E. Griffin to publish on GitHub; no named standard software licence assigned |
| RANDRAW | Alex Bar-Guy; source credits also retain Alexander Podgaetsky | Unmodified dependency used by `NG.m` | [MATLAB Central File Exchange 7309](https://uk.mathworks.com/matlabcentral/fileexchange/7309-randraw), release 1.4.0.0 (6 March 2013) | BSD 2-Clause; [full licence](sign-constrained-shrinkage-priors/normal-gamma/private/LICENSE-randraw.txt) |

## References

- Kalli, M. and Griffin, J. E. (2014). Time-varying sparsity in dynamic regression models. *Journal of Econometrics*, 178(2), 779–793. <https://doi.org/10.1016/j.jeconom.2013.10.012>
- Kowal, D. R., Matteson, D. S. and Ruppert, D. (2019). Dynamic shrinkage processes. *Journal of the Royal Statistical Society: Series B*, 81(4), 781–804. <https://doi.org/10.1111/rssb.12325>
- Griffin, J. E. and Brown, P. J. (2010). Inference with normal-gamma prior distributions in regression problems. *Bayesian Analysis*, 5(1), 171–188. <https://doi.org/10.1214/10-BA507>

Additional package and method citations are retained in the relevant component folders.
