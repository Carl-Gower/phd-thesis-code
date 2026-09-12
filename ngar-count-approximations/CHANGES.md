# Differences from the original NGAR code

This file compares the uploaded `src/NGAR.m` and `src/kf_NGAR.m` with the MATLAB implementation supplied by Kalli and Griffin.

## Count approximations

- Added optional normal and Poisson approximations for latent binomial and Poisson draws that become impractical at large counts or means.
- Retained `binornd` and `poissrnd` as the default and whenever no approximation condition applies.
- Added `n_condition` as the large-count threshold and `p_condition` to control the central and tail regions used by the binomial approximations.
- Handled the binomial endpoints explicitly and kept every approximate count within the valid count range.
- Recorded how often each ordinary, approximate or deterministic count branch is used.

## Interface and output

- Kept `NGAR` as the public entry point for both the ordinary and approximate routes.
- Accepted either row- or column-oriented responses and added checks for inputs and iteration controls.
- Added `PrintEvery` and `SaveParameters` options. These affect reporting and storage only; all parameters continue to be sampled.
- Added run settings, random-number states, elapsed time and count-branch use to `output.diagnostics`.

## Sampling implementation

- Auxiliary-count random-walk proposals use symmetric positive and negative moves.
- The persistence-parameter calculations include the beta-prior terms and logit Jacobians used by the uploaded implementation.
- Proposed latent scale paths include the required shape and rate transformations, including the equal-shape case.
- The empirical proposal covariance uses the available adaptation sample count.
- Proposed scale paths are checked before the likelihood-based acceptance step.
- Backward-sampling covariance calculations are diagonally scaled before inversion and then transformed back to the original units.
- If an empirical proposal covariance cannot be factorised, the sampler uses its independent proposal.

## Practical changes

- Identified MATLAB errors are returned for non-finite count or coefficient draws and non-positive-definite state covariances.
- Initial `lambda` and `mu` values respect the support bounds used during sampling.
- Scalar Kalman-filter divisions are evaluated directly, and the helper function name matches `kf_NGAR.m`.
- The function does not change MATLAB's global warning state.
- `lambda` and `lambdasigma` are bounded to 0.1–100; `mu` and `musigma` are bounded to 0.01–100.

## Conventions retained from the original implementation

- Persistence parameters are capped at 0.9999.
- The first coefficient keeps its special prior.
- The implemented hyperprior is proportional to `(mustar + B)^(-3)`. The paper prints `(mustar + 2B)^(-3)` while describing `B` as its mean; the code convention is retained here.
- Adaptation continues after burn-in. The tuning rules have not been redesigned for this release.

## Tests

The MATLAB tests use synthetic data to check the public interface, fixed-seed reproducibility, selected numerical identities, storage options and execution of the approximation branches. They are software checks, not a convergence study or an assessment of approximation accuracy for a particular application.
