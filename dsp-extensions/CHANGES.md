# Differences from `dsp` 0.1.0

This file describes what is different in the version uploaded here. The dynamic shrinkage process and the original Bayesian trend-filtering implementation remain the work of Kowal, Matteson and Ruppert and Daniel R. Kowal, respectively.

## State-evolution options

- Added `AR` and `DAR` coefficient-evolution options alongside the original first- and second-difference cases (`D1` and `D2`).
- Added a Gibbs update for the AR/DAR autocorrelation parameter `rho`, using either one value shared by all coefficient paths or a separate value for each predictor.
- Added the `rho_structure` argument to the regression wrapper.
- Used one explicit state operator for the evolution residuals and Gaussian precision matrices, including the stated DAR boundary convention.

## Sampling and reported quantities

- Retained both the joint coefficient draw and the path-by-path backfitting update in the same public script.
- When backfitting is selected, it is also used to initialise the coefficient paths, avoiding an initial joint factorisation.
- The univariate horseshoe observation-noise update uses the same hierarchy as the regression wrapper.
- Log-likelihood and conditional DIC calculations use the observed responses only; missing responses are sampled but are not scored as if they were data.
- Diagonal equilibration is applied before sparse Gaussian draws and reversed when the draw is returned. This changes the numerical units of the calculation, not the model.

## Interface and presentation

- Combined the sampler and helper functions required for this work into one source file.
- Added checks for dimensions, variances, iteration controls, output fields and the minimum series lengths required by each state equation.
- Added run metadata recording the resolved settings, package versions and incoming random-number state.
- Removed research-machine paths, global diagnostic assignments, unused imports and inactive simulation code from the working scripts.
- Renamed the local time-dimension argument from `T` to `n` where it improves readability. This is a change within this script, not a change to the `spam` package.

The main entry points remain `btf_reg()` and `btf()`.
