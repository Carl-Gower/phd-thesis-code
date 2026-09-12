# NGAR with optional count approximations

This folder contains my adaptation of the MATLAB implementation of the normal-gamma autoregressive (NGAR) dynamic regression model described by Kalli and Griffin (2014). I added optional approximations for large binomial and Poisson draws. The ordinary MATLAB count samplers remain the default.

> **Permission status:** Maria Kalli and Jim E. Griffin gave Carl Gower
> written permission to put the modified code on GitHub. No named standard
> software licence was assigned; see [`LICENSE.md`](LICENSE.md).

## What this version adds

- optional normal and Poisson approximations for cases in which the ordinary large-count draws become computationally impractical;
- one entry point for the ordinary and approximate sampling routes;
- selective storage of posterior histories and run diagnostics; and
- input checks and numerical safeguards listed in [`CHANGES.md`](CHANGES.md).

## Requirements

- MATLAB R2024b or a compatible release
- Statistics and Machine Learning Toolbox

## Quick start

Set MATLAB's Current Folder to `ngar-count-approximations`, then run the
self-contained example:

```matlab
run('examples/example_ngar.m')
```

It creates synthetic data and plots the posterior mean coefficient paths.
The short run demonstrates the interface; it does not establish convergence.

For your own predictor matrix `X` and response vector `y`:

```matlab
addpath('src')

rng(42)
out = NGAR(X, y, 2, 0.1, 1000, 1000, 1);
```

The full function call is:

```matlab
out = NGAR(X, y, mumean1b, mulambdastar, burnin, numbofits, every, ...
    n_condition, p_condition, options);
```

The first seven inputs are required. The final thresholds and `options` are optional:

| Input | Meaning |
| --- | --- |
| `X` | `T`-by-`p` predictor matrix, including an intercept as its first column if required |
| `y` | Response vector of length `T` |
| `mumean1b` | Positive scale `B` in the implemented hyperprior proportional to `(mustar + B)^(-3)` |
| `mulambdastar` | Positive mean of the exponential prior for `lambdastar` |
| `burnin` | Number of burn-in iterations |
| `numbofits` | Number of posterior draws to save |
| `every` | Number of iterations between saved draws |
| `n_condition` | Large-count threshold at which an approximation may be used |
| `p_condition` | Probability-range control for choosing between the normal and Poisson binomial approximations; between 0 and 0.5 |
| `options` | Optional scalar structure controlling progress output and saved parameters |

`X` must be a finite `T`-by-`p` matrix and `y` a vector of length `T`. Include
an intercept explicitly as the first column of `X`; the first coefficient has
the special prior used by the original implementation.

To use the approximation rules, supply both thresholds:

```matlab
rng(42)
outApprox = NGAR(X, y, 2, 0.1, 1000, 1000, 1, 100, 0.25);
```

When a relevant count or Poisson mean is below `n_condition`, the ordinary MATLAB draw is retained. At or above the threshold, `p_condition` helps select a normal approximation near the centre of a binomial distribution or a Poisson approximation in its tails. Omitting both thresholds uses `binornd` and `poissrnd` throughout. The total number of sampler iterations is `burnin + every*numbofits`.

## Output and run options

The default output contains 12 parameter histories and a `diagnostics` struct.
With `S = numbofits`, the history dimensions are:

| Fields | Dimensions |
| --- | --- |
| `beta`, `Psi` | `T`-by-`p`-by-`S` |
| `sigmasq` | `T`-by-`S` |
| `lambda`, `mu`, `rhobeta`, `rho` | `p`-by-`S` |
| `lambdasigma`, `musigma`, `rhosigma`, `lambdastar`, `mustar` | `1`-by-`S` |

For example, `mean(out.beta, 3)` gives the posterior mean coefficient paths.

An optional final struct controls progress printing and output storage:

```matlab
options.PrintEvery = 0;
options.SaveParameters = {'beta', 'Psi'};

out = NGAR(X, y, 2, 0.1, 1000, 1000, 1, options);
```

All parameters are still sampled when only selected histories are saved. See
[`docs/DISPLAY_DIAGNOSTICS.md`](docs/DISPLAY_DIAGNOSTICS.md) for the complete
option and diagnostic definitions.

## Tests

Run the regression tests from this folder with:

```matlab
results = runtests('tests');
assertSuccess(results)
```

The supplied tests cover the interface, reproducibility, numerical identities, storage options and execution of the approximation branches.

## Important limitations

- Approximation mode uses approximate proposal draws. Its invariant distribution has not been established as identical to that of the ordinary route, so it should be treated as approximate MCMC.
- The regression tests are software checks. They do not establish convergence,
  posterior accuracy or an appropriate approximation threshold for a new data
  set.
- The parameter bounds and adaptation scheme inherited from the research
  implementation are documented in [`CHANGES.md`](CHANGES.md).

## Provenance and citation

The relationship to the original research code is set out in [`PROVENANCE.md`](PROVENANCE.md). [`CHANGES.md`](CHANGES.md) gives a direct comparison between the uploaded version and the original implementation.

Kalli, M. and Griffin, J. E. (2014). Time-varying sparsity in dynamic
regression models. *Journal of Econometrics*, 178, 779–793.
[https://doi.org/10.1016/j.jeconom.2013.10.012](https://doi.org/10.1016/j.jeconom.2013.10.012)

If you use the optional count approximations or the uploaded implementation, please cite this repository as well as Kalli and Griffin (2014). GitHub displays the repository citation from the root [`CITATION.cff`](../CITATION.cff).

## Licence

Published on GitHub with written permission from Maria Kalli and Jim E.
Griffin. No named standard software licence is asserted for this component;
see [`LICENSE.md`](LICENSE.md).
