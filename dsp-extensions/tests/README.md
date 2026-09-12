# Validation suite

`validate_dsp.R` is a reproducible synthetic test suite for `R/DSP_rho.R`. It uses base R assertions, so it does not add a test-framework dependency.

From the component root, run:

```sh
Rscript tests/validate_dsp.R
```

To check another copy of the source file:

```sh
Rscript tests/validate_dsp.R path/to/DSP_rho.R
```

The suite covers:

- independently constructed state operators and dense posterior precisions;
- conditional Gaussian means and covariances;
- simulated moments from the sparse joint draw;
- invariance of one random-order backfitting sweep;
- AR and DAR draws for `rho`, including tail cases;
- short complete MCMC runs across all evolution equations, shrinkage priors, observation-variance choices, coefficient samplers and `rho` structures;
- missing responses, collinear and zero predictors, minimum supported lengths and one-predictor models;
- invalid inputs, output selection and seeded replay.

The Monte Carlo checks use fixed seeds, but results can differ if numerical libraries or dependency implementations change. A deliberately very short test series may make the base R `arima()` initialisation issue a warning; that initialisation is separate from posterior sampling of `rho`.

These tests address implementation consistency and numerical behaviour. They do not assess posterior convergence for a substantive application or reproduce the thesis analysis.
