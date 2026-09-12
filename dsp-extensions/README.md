# Dynamic shrinkage process extensions

This folder contains the standalone R implementation I used to extend the Bayesian trend-filtering samplers in Daniel Kowal's [`dsp`](https://github.com/drkowal/dsp) package. The underlying dynamic shrinkage process is from Kowal, Matteson and Ruppert (2019). My additions concern the evolution equation for the time-varying coefficients and the way its autocorrelation parameter is handled.

The main script supports four evolution specifications (`AR`, `D1`, `DAR` and `D2`). For the two specifications involving an autocorrelation parameter, `rho` can be shared by all coefficient paths or estimated separately for each predictor. Both a joint coefficient draw and a backfitting update are available.

This is research code rather than an R package, and it is not an official release of `dsp`.

## Quick start

The code was validated with R 4.5.1 and these package versions:

| Package | Version used for validation |
| --- | ---: |
| `Matrix` | 1.7-3 |
| `spam` | 2.11-1 |
| `stochvol` | 3.2.5 |
| `BayesLogit` | 2.1 |
| `truncdist` | 1.0-2 |

Install the dependencies once:

```r
install.packages(c("Matrix", "spam", "stochvol", "BayesLogit", "truncdist"))
```

Then, from this folder, run the synthetic example:

```sh
Rscript examples/example_dsp.R
```

The example deliberately uses a short chain so that it finishes quickly. It demonstrates the interface; it is not a recommendation for empirical work.

To use the sampler in another script:

```r
source("R/DSP_rho.R")

set.seed(1234)
fit <- btf_reg(
  y,
  X,
  evol_error = "DHS",
  differencing_option = "DAR",
  rho_structure = "variable_specific",
  use_backfitting = TRUE,
  nsave = 1000,
  nburn = 1000,
  nskip = 4
)
```

`X` must be a finite numeric matrix with one row per observation. Include a column of ones if an intercept is required. Predictors are not scaled automatically. `y` must be a numeric vector with at least two observed values and positive sample variance. The sampler needs at least four time points, or five for `D2`.

## Evolution specifications

For coefficient path `j`, the state residuals are:

| Option | State residual |
| --- | --- |
| `AR` | `beta[t,j] - rho[j] * beta[t-1,j]` |
| `D1` | `beta[t,j] - beta[t-1,j]` |
| `DAR` | `beta[t,j] - (1 + rho[j]) * beta[t-1,j] + rho[j] * beta[t-2,j]` |
| `D2` | `beta[t,j] - 2 * beta[t-1,j] + beta[t-2,j]` |

The first DAR residual is defined as `beta[2,j] - (1-rho[j]) * beta[1,j]`. This boundary equation is a modelling convention used consistently by both coefficient samplers. D2 gives the first two states separate proper priors; the other specifications give the first state a proper prior.

For `AR` and `DAR`, the public wrappers use a `Uniform(0, 0.99)` prior for `rho`. Set `rho_structure = "shared"` for one value across all predictors, or `"variable_specific"` for one value per path. D1 and D2 do not use `rho`.

The available evolution-error priors are dynamic horseshoe (`DHS`), horseshoe (`HS`), Bayesian lasso (`BL`), stochastic volatility (`SV`) and normal-inverse-gamma (`NIG`). The default is `DHS`.

## Main functions

- `btf_reg()` fits the time-varying regression model.
- `btf()` fits the univariate smoothing model.
- `sampleBTF_reg()` draws all coefficient paths jointly.
- `sampleBTF_reg_backfit()` updates the paths one at a time.
- `sample_AR1_param()` draws a shared or predictor-specific `rho` conditional on the coefficient paths.

The fitted object contains the requested draws and always includes `loglike`. A `run_info` attribute records the resolved settings, R version, dependency versions and the incoming random-number state when one exists.

For a regression with `n` observations and `p` predictors, `fit$beta` has
dimensions `nsave`-by-`n`-by-`p`. For example,
`apply(fit$beta, c(2, 3), median)` gives an `n`-by-`p` matrix of posterior
median coefficient paths. The `mu`, `yhat` and `obs_sigma_t2` histories each
have dimensions `nsave`-by-`n`.

## Missing values and one-step forecasts

Missing values in `y` are sampled within the model. An internal `NA` is therefore treated as a missing observation. To obtain posterior predictive draws for the next response, append the next predictor row to `X`, append `NA` to `y`, and request `yhat`:

```r
y_forecast <- c(y, NA_real_)
X_forecast <- rbind(X, X_next)

fit <- btf_reg(
  y_forecast,
  X_forecast,
  mcmc_params = c("beta", "yhat"),
  nsave = 1000,
  nburn = 1000
)

next_y_draws <- fit$yhat[, length(y_forecast)]
```

`X_next` must have the same columns, in the same order, as `X`. Appending several rows with `NA` gives multi-step predictive draws, conditional on the supplied future predictor rows.

## Validation

Run the full synthetic validation suite from this folder with:

```sh
Rscript tests/validate_dsp.R
```

The suite checks the state operators against independent dense equations, the means and covariances of the Gaussian draws, a backfitting sweep, the conditional distribution of `rho`, input handling, reproducibility under a fixed seed and short end-to-end fits across the supported options. See [`tests/README.md`](tests/README.md) for scope and caveats.

## Points to bear in mind

- With a fully observed `y`, the returned coefficient paths are smoothing draws conditional on the full supplied series. Forecasts require one or more appended `NA` responses and the corresponding future predictor rows, as shown above.
- The example and validation suite use synthetic data. They do not reproduce the thesis results or establish MCMC convergence for a new application.
- Backfitting reduces the size of each Cholesky factorisation, but it may mix more slowly when predictors are strongly correlated.
- Extremely disparate variance scales can still make a precision matrix numerically difficult. The sampler reports a failed factorisation rather than adding a small ridge term.
- The source is kept as one script to preserve the workflow used in the research. Sourcing it adds its functions to the current R environment.

## Original work and citation

The method and original software should be credited to:

> Kowal, D. R., Matteson, D. S. and Ruppert, D. (2019). Dynamic shrinkage processes. *Journal of the Royal Statistical Society: Series B*, 81(4), 781-804. <https://doi.org/10.1111/rssb.12325>

The exact upstream commit used for comparison was [`d5ea2548cb251c078b79f57a5a9e273d2650f4bd`](https://github.com/drkowal/dsp/tree/d5ea2548cb251c078b79f57a5a9e273d2650f4bd). [`PROVENANCE.md`](PROVENANCE.md) records the source correspondence, and [`CHANGES.md`](CHANGES.md) distinguishes my changes from the upstream implementation.

If you use the AR/DAR specifications, the `rho` extensions or this standalone implementation, please also cite this repository. GitHub displays its citation from the root [`CITATION.cff`](../CITATION.cff).

## License

Because this is a modified version of GPL-2 code, this component is distributed under the [GNU General Public License version 2 only](LICENSE). The licence applies to the modified work in this folder. The R packages on which it depends retain their own licences.
