# Sign-constrained static shrinkage priors

This folder collects the changes I made to several existing Bayesian regression samplers so that selected slopes can be restricted to positive or negative values. The motivation was to impose the economic sign restrictions used by Campbell and Thompson (2008), but the interface is not tied to that application.

These are source overlays for specific upstream versions, not new packages. Keeping only the changed files makes my contribution clearer than republishing several complete packages. The supplied helper applies an overlay to an upstream source folder so that the files do not have to be replaced individually.

## Implementations

| Folder | Upstream version | Constraint arguments | Index convention |
| --- | --- | --- | --- |
| [`bayesreg/`](bayesreg/) | `bayesreg` 1.3 | `pos_con`, `neg_con` | Columns of the predictor matrix after the intercept is removed |
| [`r2d2/`](r2d2/) | R2D2 commit `e7346399` | `pos_con`, `neg_con` | Columns of `x` |
| [`boom-spike-slab/`](boom-spike-slab/) | `Boom` 0.9.16 and `BoomSpikeSlab` 1.2.7 | `pos_con`, `neg_con` | Columns of the full model matrix, including its intercept column |
| [`normal-gamma/`](normal-gamma/) | Code supplied by Jim E. Griffin | `pos_con`, `neg_con` | Columns of `data`; the separately sampled intercept is not counted |

The two index vectors must contain distinct one-based integer positions. An empty vector leaves that side unconstrained. Because the coefficient distributions are continuous, a zero boundary corresponds to positive or negative draws with probability one, apart from numerical boundary cases and excluded spike-and-slab coefficients, which remain zero.

## Installing the R overlays

Use R 4.5.0 or newer and a working source-package toolchain. `Boom` and `BoomSpikeSlab` compile C++17 code, so the first installation can take several minutes. On macOS, use the [R for macOS toolchain instructions](https://mac.r-project.org/tools/); on Windows, install the [Rtools version matching your R release](https://cran.r-project.org/bin/windows/Rtools/). The workflow has been tested on macOS with R 4.5.1; Windows and Linux have not been tested here.

Set R's working directory to this `sign-constrained-shrinkage-priors` folder, then run:

```r
source("install_overlays.R")
source("use_overlays.R")
source("examples/example_constraints.R")
```

The installer downloads the exact upstream sources, checks their checksums, applies the modified files, installs missing dependencies, and runs the installed-package checks. It builds `Boom` before `BoomSpikeSlab` and makes sure the latter uses the modified `Boom` headers and library.

The modified packages are installed in `.r-library/` inside this folder, separately for each R series and platform. Your usual package library is left intact. Downloads, build files and installation logs go in `.r-build/`; both folders are excluded from Git. Dependency versions are taken from your existing libraries or CRAN, so this is not a complete lockfile for the R environment.

In each new R session, run `source("use_overlays.R")` **before loading any of these packages**. If an ordinary version is already loaded, the helper will explain that you need to restart R once and run it again. There is no need to delete or uninstall your existing packages. Changing the library search path cannot replace a namespace or compiled library already loaded in R.

You can also install and run the examples in fresh sessions from a terminal:

```sh
Rscript --vanilla install_overlays.R
Rscript --vanilla examples/example_constraints.R
```

If installation fails, the installer prints the end of the error log and its full location. A missing compiler must be installed before trying again. Running only `apply_overlay.R` copies source files; it does not install or activate a package.

### Applying an individual overlay manually

Use this route only if you want to manage the source folders and installation yourself. Obtain `bayesreg` 1.3, `Boom` 0.9.16, `BoomSpikeSlab` 1.2.7 or the recorded [R2D2 commit](https://github.com/yandorazhang/R2D2/tree/e734639929abb60e616c114ac7fe4e2beb5c7f9d). The exact download URLs, including CRAN archive fallbacks, are recorded in [`install_overlays.R`](install_overlays.R).

For example, apply the `bayesreg` changes to a freshly unpacked source folder:

```sh
Rscript apply_overlay.R bayesreg /path/to/bayesreg
```

The component names are `bayesreg`, `r2d2`, `boom` and `boomspikeslab`. The helper checks the package name, version and source files before copying, rejects installed package folders, and labels modified versions with a `.9000` suffix. Reapplying it to the same source folder is supported. You must then install the source package and its dependencies yourself. When building `BoomSpikeSlab`, both R's library path and the `R_LIBS` environment variable must select the modified `Boom` first; the automatic installer handles this.

## Minimal calls

```r
# bayesreg: x1 positive and x2 negative (the intercept is not counted)
fit_br <- bayesreg::bayesreg(
  y ~ x1 + x2, data = dat, prior = "horseshoe",
  n.samples = 1000, burnin = 1000, thin = 1, n.cores = 1,
  pos_con = 1, neg_con = 2
)

# R2D2: columns 1 and 2 of X
fit_dl <- R2D2::dl(X, y, mcmc.n = 2000, print = FALSE,
                   pos_con = 1, neg_con = 2)

# BoomSpikeSlab: columns of model.matrix(y ~ x1 + x2)
fit_ss <- BoomSpikeSlab::lm.spike(
  y ~ x1 + x2, niter = 2000, data = dat,
  model.options = BoomSpikeSlab::SsvsOptions(adaptive.cutoff = Inf),
  pos_con = 2, neg_con = 3
)
```

Runnable synthetic examples are in [`examples/`](examples/).

## Supported paths

- The `bayesreg` constraints are only available for the Gaussian, Laplace, Student-t and logistic scale-mixture samplers.
- The `BoomSpikeSlab` constraints are implemented only for Gaussian regression with the classic SSVS sampler. The public R wrapper stops if constraints are requested with adaptive SSVS, ODA or Student errors.
- The `Boom` coefficient update uses a coordinate-wise Gibbs sweep for constrained coefficients in the supported classic SSVS route.
- The Normal-Gamma MATLAB implementation requires Statistics and Machine Learning Toolbox. Its `randraw.m` dependency is included in `normal-gamma/private/` with its BSD licence; no separate download is needed.

## Scope of the tests

- The supplied checks confirm that the supported routes run and that saved coefficient draws have the requested signs on small synthetic examples.
- They do not establish posterior equivalence with a separate implementation or assess MCMC convergence and mixing for a substantive application.
- As with any MCMC analysis, convergence and effective sample sizes should be assessed for the data and settings actually used.

## Tests

After installing the overlays, run:

```sh
Rscript --vanilla tests/smoke_installed.R
```

The test checks argument validation and sign compliance on small synthetic problems. See [`tests/README.md`](tests/README.md) for the exact scope.

## Provenance, licences and permission

[`PROVENANCE.md`](PROVENANCE.md) records the upstream versions and changed files. Each R component retains its upstream licence. The Normal-Gamma source is published on GitHub with written permission from Jim E. Griffin. No named standard software licence was assigned; see [`normal-gamma/LICENSE.md`](normal-gamma/LICENSE.md).

If you use the sign-constraint modifications, please cite this repository as well as the relevant original method and software. GitHub displays the repository citation from the root [`CITATION.cff`](../CITATION.cff).

## References

- Campbell, J. Y. and Thompson, S. B. (2008). Predicting excess stock returns out of sample: can anything beat the historical average? *Review of Financial Studies*, 21(4), 1509–1531. <https://doi.org/10.1093/rfs/hhm055>
- Griffin, J. E. and Brown, P. J. (2010). Inference with normal-gamma prior distributions in regression problems. *Bayesian Analysis*, 5(1), 171–188. <https://doi.org/10.1214/10-BA507>

The package-specific papers and authors remain documented in the upstream packages.
