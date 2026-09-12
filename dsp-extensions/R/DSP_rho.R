# Bayesian trend filtering with optional AR/DAR state evolution
#
# Standalone adaptation of dsp 0.1.0 by Daniel R. Kowal.
# Based on Kowal, Matteson and Ruppert (2019), Dynamic Shrinkage Processes,
# https://doi.org/10.1111/rssb.12325
# Upstream source: https://github.com/drkowal/dsp
# Reference commit: d5ea2548cb251c078b79f57a5a9e273d2650f4bd
#
# PhD extensions by Carl Gower, 2025:
#   * added AR and DAR state-evolution options;
#   * added shared and predictor-specific autocorrelation parameters;
#   * combined the regression and component samplers in one script; and
#   * aligned the state equations across the joint and backfitting samplers.
# A fuller record of the changes is in ../CHANGES.md and ../PROVENANCE.md.
#
# This modified work is free software distributed under version 2 only of the
# GNU General Public License. It is provided without warranty. See ../LICENSE.
# SPDX-License-Identifier: GPL-2.0-only
#
# Model conventions:
#   n = number of time points, p = number of predictors.
#   X must contain any desired intercept. No automatic scaling is performed.
#   The DAR boundary equation is documented above .state_operator().
#
# See ../README.md and ../examples/example_dsp.R for a reproducible example.

# Check all dependencies at once; sourcing never installs packages or runs a fit.
local({
  required <- c("stochvol", "Matrix", "spam", "BayesLogit", "truncdist")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages before sourcing: ", paste(missing, collapse = ", "))
})
suppressPackageStartupMessages({
  library(stochvol)
  library(Matrix)
  library(spam)
  library(BayesLogit)
  library(truncdist)
})

# ============================================================================
# User-facing samplers and AR/DAR extensions
# ============================================================================

#' MCMC sampler for Bayesian trend-filtering regression
#'
#' Extends `dsp::btf_reg()` with AR and DAR state evolution. For AR and DAR,
#' `rho_structure` selects either one shared autocorrelation parameter or one
#' parameter per predictor. D1 and D2 do not use `rho`.
#'
#' @param y Numeric length-n response vector; NA values are imputed. At least
#'   two observed values with positive variance and four time points (five for
#'   D2) are required by the retained evolution-prior initialisation/updates.
#' @param X Finite numeric n by p predictor matrix, or NULL for the univariate
#'   model. Include an intercept explicitly; predictors are not rescaled.
#' @param evol_error Evolution-error prior: DHS, HS, BL, SV, or NIG.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho_structure Either shared or variable_specific.
#' @param useObsSV Whether to use stochastic volatility for observation errors.
#' @param nsave Positive integer number of draws to save; at least two for DIC.
#' @param nburn Number of burn-in draws.
#' @param nskip Number of draws skipped between saved draws.
#' @param mcmc_params Character vector or list of names: mu, yhat, beta,
#'   evol_sigma_t2, obs_sigma_t2, dhs_phi, dhs_mean, or rho.
#' @param use_backfitting Whether to use conditional path sweeps, including
#'   initialisation. FALSE draws all paths jointly; TRUE uses smaller Cholesky
#'   solves and can mix more slowly when predictors are strongly correlated.
#' @param computeDIC Whether to calculate conditional DIC on observed responses.
#' @param verbose Whether to report timing information.
#'
#' @return Named list of requested draws: beta and evol_sigma_t2 are
#'   nsave by n by p; mu, yhat and obs_sigma_t2 are nsave by n. DHS parameters
#'   are nsave by p; rho is length nsave (shared) or nsave by p. Inapplicable
#'   requested parameters remain NULL. loglike is always returned. DIC and p_d
#'   each contain two versions when requested; run_info is a reproducibility
#'   attribute. X = NULL returns the univariate btf() output.
#' @note AR/DAR use a Uniform(0, 0.99) rho prior. The returned paths are
#'   smoothing draws conditional on the whole supplied series.
btf_reg <- function(
    y,
    X = NULL,
    evol_error = "DHS",
    differencing_option = "D1",
    rho_structure = c("shared", "variable_specific"),
    useObsSV = FALSE,
    nsave = 1000,
    nburn = 1000,
    nskip = 4,
    mcmc_params = list(
      "mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2",
      "dhs_phi", "dhs_mean", "rho"),
    use_backfitting = FALSE,
    computeDIC = TRUE,
    verbose = TRUE) {

  # Validate controls before allocating arrays or entering the MCMC.
  call <- match.call()
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  evol_error <- match.arg(toupper(evol_error), c("DHS", "HS", "BL", "SV", "NIG"))
  differencing_option <- match.arg(differencing_option, c("AR", "D1", "DAR", "D2"))
  mcmc_params <- .check_mcmc(y, differencing_option,
    list(nsave = nsave, nburn = nburn, nskip = nskip, useObsSV = useObsSV,
         computeDIC = computeDIC, verbose = verbose, use_backfitting = use_backfitting),
    mcmc_params, regression = TRUE)
  rho_structure <- match.arg(rho_structure)

  if (is.null(X)) {
    # With no design matrix, use the corresponding univariate smoothing model.
    mcmc_params <- setdiff(mcmc_params, "beta")
    return(btf(
      y = y,
      evol_error = evol_error,
      differencing_option = differencing_option,
      useObsSV = useObsSV,
      nsave = nsave,
      nburn = nburn,
      nskip = nskip,
      mcmc_params = mcmc_params,
      computeDIC = computeDIC,
      verbose = verbose
    ))
  }

  .check_data(y, X, allow_missing = TRUE)
  n <- length(y)
  # Missing responses are initialised once, then imputed within every sweep.
  observed <- !is.na(y)
  is.missing <- which(is.na(y))
  any.missing <- length(is.missing) > 0

  if (any.missing) y[is.missing] <- mean(y, na.rm = TRUE)

  # Backfitting only needs per-predictor likelihood terms.
  XtX <- if (use_backfitting) NULL else build_XtX(X)
  p <- ncol(X)

  sigma_e <- sd(y, na.rm = TRUE)
  sigma_et <- rep(sigma_e, n)

  if (differencing_option %in% c("AR", "DAR")) {
    if (rho_structure == "shared") {
      rho <- sample_AR1_param(
        beta = NULL,
        differencing_option = differencing_option,
        prior_type = "uniform"
      )
    } else {
      rho <- vapply(seq_len(p), function(j) {
        sample_AR1_param(
          beta = NULL,
          differencing_option = differencing_option,
          prior_type = "uniform"
        )
      }, numeric(1))
    }
  } else {
    rho <- NULL
  }

  # Initialisation must respect the selected coefficient sampler too.
  # Starting from zero is safe here: one proper conditional draw produces beta.
  initial_variances <- matrix(0.01 * sigma_et^2, nrow = n, ncol = p)
  chol0 <- NULL
  if (use_backfitting) {
    beta <- sampleBTF_reg_backfit(y, X, matrix(0, n, p), sigma_et^2,
                                  initial_variances, differencing_option, rho)
  } else {
    chol0 <- initCholReg.spam(sigma_et^2, initial_variances, XtX,
                             differencing_option, rho)
    beta <- sampleBTF_reg(y, X, sigma_et^2, initial_variances, XtX,
                          differencing_option, rho, chol0)
  }

  mu <- rowSums(X * beta)
  omega <- compute_omega(beta, differencing_option, rho)
  beta0 <- compute_beta0(beta, differencing_option)

  evolParams <- initEvolParams(omega, evol_error = evol_error)
  evolParams0 <- initEvol0(beta0, commonSD = FALSE)

  if (useObsSV) {
    svParams <- initSV(y - mu)
    sigma_et <- as.numeric(svParams$sigma_wt)
  }

  # Allocate requested draws; DIC additionally needs mean and variance draws.
  mcmc_output <- vector("list", length(mcmc_params))
  names(mcmc_output) <- mcmc_params

  if (!is.na(match("mu", mcmc_params)) || computeDIC) {
    post_mu <- array(NA, c(nsave, n))
  }
  if (!is.na(match("yhat", mcmc_params))) {
    post_yhat <- array(NA, c(nsave, n))
  }
  if (!is.na(match("beta", mcmc_params))) {
    post_beta <- array(NA, c(nsave, n, p))
  }
  if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) {
    post_obs_sigma_t2 <- array(NA, c(nsave, n))
  }
  if (!is.na(match("evol_sigma_t2", mcmc_params))) {
    post_evol_sigma_t2 <- array(NA, c(nsave, n, p))
  }
  if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == "DHS") {
    post_dhs_phi <- array(NA, c(nsave, p))
  }
  if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == "DHS") {
    post_dhs_mean <- array(NA, c(nsave, p))
  }
  if (differencing_option %in% c("AR", "DAR")) {
    if (rho_structure == "shared") {
      post_rho <- numeric(nsave)
    } else {
      post_rho <- matrix(NA_real_, nrow = nsave, ncol = p)
    }
  }
  post_loglike <- numeric(nsave)

  nstot <- nburn + (nskip + 1) * nsave
  skipcount <- 0
  isave <- 0

  if (verbose) timer0 <- proc.time()[3]

  # One Gibbs sweep: missing y, rho, coefficient paths, shrinkage, then noise.
  for (nsi in seq_len(nstot)) {
    if (any.missing) {
      y[is.missing] <- mu[is.missing] +
        sigma_et[is.missing] * rnorm(length(is.missing))
    }

    # Initial-state priors precede the evolution-error variances in H beta.
    initial_rows <- if (differencing_option == "D2") 2L else 1L
    evol_sigma_t2 <- rbind(
      matrix(evolParams0$sigma_w0^2, nrow = initial_rows),
      evolParams$sigma_wt^2
    )

    if (differencing_option %in% c("AR", "DAR")) {
      if (rho_structure == "shared") {
        rho <- sample_AR1_param(
          beta = beta,
          differencing_option = differencing_option,
          evol_sigma_t2 = evol_sigma_t2,
          rho_current = rho,
          prior_type = "uniform"
        )
      } else {
        for (j in seq_len(p)) {
          rho[j] <- sample_AR1_param(
            beta = beta[, j, drop = FALSE],
            differencing_option = differencing_option,
            evol_sigma_t2 = evol_sigma_t2[, j, drop = FALSE],
            rho_current = rho[j],
            prior_type = "uniform"
          )
        }
      }
    }

    # Choose one joint draw or a sweep of conditional path draws.
    if (use_backfitting) {
      beta <- sampleBTF_reg_backfit(
        y,
        X,
        beta,
        obs_sigma_t2 = sigma_et^2,
        evol_sigma_t2 = evol_sigma_t2,
        differencing_option = differencing_option,
        rho = rho
      )
    } else {
      beta <- sampleBTF_reg(
        y,
        X,
        obs_sigma_t2 = sigma_et^2,
        evol_sigma_t2 = evol_sigma_t2,
        XtX = XtX,
        differencing_option = differencing_option,
        rho = rho,
        chol0 = chol0
      )
    }

    # Recompute residuals under the current rho before updating shrinkage.
    mu <- rowSums(X * beta)
    omega <- compute_omega(beta, differencing_option, rho)
    beta0 <- compute_beta0(beta, differencing_option)

    evolParams0 <- sampleEvol0(
      beta0,
      evolParams0,
      A = 1,
      commonSD = FALSE
    )

    # Update shrinkage scales and either time-varying or constant noise.
    if (useObsSV) {
      evolParams <- sampleEvolParams(
        omega,
        evolParams,
        1 / sqrt(n * p),
        evol_error
      )
      svParams <- sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et <- as.numeric(svParams$sigma_wt)
    } else {
      evolParams <- sampleEvolParams(
        omega,
        evolParams,
        sigma_e / sqrt(n * p),
        evol_error
      )

      if (evol_error == "DHS") {
        sigma_e <- uni.slice(sigma_e, g = function(x) {
          -(n + 2) * log(x) -
            0.5 * sum((y - mu)^2, na.rm = TRUE) / x^2 -
            log(1 + (sqrt(n * p) * exp(evolParams$dhs_mean0 / 2) / x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if (evol_error == "HS") {
        sigma_e <- 1 / sqrt(rgamma(
          n = 1,
          shape = n / 2,
          rate = sum((y - mu)^2, na.rm = TRUE) / 2
        ))
      }
      if (evol_error == "BL") {
        sigma_e <- 1 / sqrt(rgamma(
          n = 1,
          shape = n / 2 + length(evolParams$tau_j) / 2,
          rate = sum((y - mu)^2, na.rm = TRUE) / 2 +
            n * p * sum((omega / evolParams$tau_j)^2) / 2
        ))
      }
      if (evol_error %in% c("NIG", "SV")) {
        sigma_e <- 1 / sqrt(rgamma(
          n = 1,
          shape = n / 2,
          rate = sum((y - mu)^2, na.rm = TRUE) / 2
        ))
      }

      sigma_et <- rep(sigma_e, n)
    }

    # Retain every (nskip + 1)-th draw after burn-in.
    if (nsi > nburn) {
      skipcount <- skipcount + 1

      if (skipcount > nskip) {
        isave <- isave + 1

        if (!is.na(match("mu", mcmc_params)) || computeDIC) {
          post_mu[isave, ] <- mu
        }
        if (!is.na(match("yhat", mcmc_params))) {
          post_yhat[isave, ] <- mu + sigma_et * rnorm(n)
        }
        if (!is.na(match("beta", mcmc_params))) {
          post_beta[isave, , ] <- beta
        }
        if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) {
          post_obs_sigma_t2[isave, ] <- sigma_et^2
        }
        if (!is.na(match("evol_sigma_t2", mcmc_params))) {
          post_evol_sigma_t2[isave, , ] <- rbind(
            matrix(evolParams0$sigma_w0^2, nrow = initial_rows),
            evolParams$sigma_wt^2
          )
        }
        if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == "DHS") {
          post_dhs_phi[isave, ] <- evolParams$dhs_phi
        }
        if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == "DHS") {
          post_dhs_mean[isave, ] <- evolParams$dhs_mean
        }
        if (differencing_option %in% c("AR", "DAR")) {
          if (rho_structure == "shared") {
            post_rho[isave] <- rho
          } else {
            post_rho[isave, ] <- rho
          }
        }

        # Score observed responses only; imputations are latent variables.
        post_loglike[isave] <- sum(dnorm(
          y[observed],
          mean = as.numeric(mu)[observed],
          sd = sigma_et[observed],
          log = TRUE
        ))
        skipcount <- 0
      }
    }

    if (verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if (!is.na(match("mu", mcmc_params))) mcmc_output$mu <- post_mu
  if (!is.na(match("yhat", mcmc_params))) mcmc_output$yhat <- post_yhat
  if (!is.na(match("beta", mcmc_params))) mcmc_output$beta <- post_beta
  if (!is.na(match("obs_sigma_t2", mcmc_params))) {
    mcmc_output$obs_sigma_t2 <- post_obs_sigma_t2
  }
  if (!is.na(match("evol_sigma_t2", mcmc_params))) {
    mcmc_output$evol_sigma_t2 <- post_evol_sigma_t2
  }
  if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == "DHS") {
    mcmc_output$dhs_phi <- post_dhs_phi
  }
  if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == "DHS") {
    mcmc_output$dhs_mean <- post_dhs_mean
  }
  if (differencing_option %in% c("AR", "DAR") &&
      !is.na(match("rho", mcmc_params))) {
    mcmc_output$rho <- post_rho
  }

  mcmc_output$loglike <- post_loglike

  if (computeDIC) {
    # Conditional DIC follows dsp's two p_d definitions, using observed y only.
    loglike_hat <- sum(dnorm(
      y[observed],
      mean = colMeans(post_mu)[observed],
      sd = colMeans(sqrt(post_obs_sigma_t2))[observed],
      log = TRUE
    ))
    p_d <- c(
      2 * (loglike_hat - mean(post_loglike)),
      2 * var(post_loglike)
    )
    mcmc_output$DIC <- -2 * loglike_hat + 2 * p_d
    mcmc_output$p_d <- p_d
  }

  if (verbose) {
    print(paste("Total time: ", round(proc.time()[3] - timer0), "seconds"))
  }

  # Record reproducibility information alongside, rather than among, draws.
  attr(mcmc_output, "run_info") <- list(
    call = call, n = n, p = p, observed_n = sum(observed),
    rho_structure = rho_structure, use_backfitting = use_backfitting,
    differencing_option = differencing_option, evol_error = evol_error,
    nsave = nsave, nburn = nburn, nskip = nskip,
    useObsSV = useObsSV, computeDIC = computeDIC, mcmc_params = mcmc_params,
    rng_kind = RNGkind(), rng_state_at_entry = rng_state,
    R_version = R.version.string,
    package_versions = vapply(c("Matrix", "spam", "stochvol", "BayesLogit", "truncdist"),
                              function(pkg) as.character(utils::packageVersion(pkg)), character(1)))
  mcmc_output
}

# ============================================================================
# State equations, autocorrelation updates, and Gaussian coefficient samplers
# ============================================================================

#' Validate a scalar-or-vector autocorrelation specification
#' @param rho Scalar or length-p numeric vector, with absolute values below one.
#' @param p Number of coefficient paths.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @return Length-p vector for AR/DAR; NULL for D1/D2.
#' @keywords internal
.rho_vector <- function(rho, p, differencing_option) {
  if (!differencing_option %in% c("AR", "DAR")) return(NULL)
  if (!is.numeric(rho) || !length(rho) %in% c(1L, p) ||
      any(!is.finite(rho)) || any(abs(rho) >= 1)) {
    stop("rho must be a finite scalar or length-p vector with abs(rho) < 1.")
  }
  rep(as.numeric(rho), length.out = p)
}

#' Validate variances before forming a Gaussian precision matrix
#' @param x Numeric variance vector or matrix.
#' @param label Argument name used in an error message.
#' @param allow_inf Whether positive infinity can represent zero precision.
#' @return Invisibly, NULL; invalid inputs raise an informative error.
#' @keywords internal
.check_variances <- function(x, label, allow_inf = FALSE) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(x <= 0) ||
      (!allow_inf && any(!is.finite(x)))) {
    stop(label, " must contain strictly positive ",
         if (allow_inf) "values (Inf is allowed)." else "finite values.")
  }
  invisible(NULL)
}

#' Validate the observed response and regression design
#' @param y Numeric response vector.
#' @param X Numeric n by p matrix, or NULL.
#' @param allow_missing Whether NA responses can be imputed by the MCMC wrapper.
#' @return Invisibly, NULL; invalid inputs raise an informative error.
#' @keywords internal
.check_data <- function(y, X = NULL, allow_missing = FALSE) {
  if (!is.numeric(y) || !is.null(dim(y)) || !length(y) ||
      any(is.infinite(y)) || any(is.nan(y)) ||
      (!allow_missing && anyNA(y))) stop("y must be a numeric vector of finite observations", if (allow_missing) " or NA." else ".")
  if (!is.null(X) && (!is.matrix(X) || !is.numeric(X) ||
      nrow(X) != length(y) || ncol(X) < 1L || any(!is.finite(X)))) {
    stop("X must be a finite numeric matrix with length(y) rows and at least one column.")
  }
  invisible(NULL)
}

#' Validate MCMC controls and requested output fields
#' @param y Response vector; NA values are allowed.
#' @param differencing_option State evolution specification.
#' @param controls Named list containing nsave, nburn, nskip and logical switches.
#' @param mcmc_params Character vector or list of output names.
#' @param regression Whether beta draws are available.
#' @return Unique character vector of requested output names.
#' @keywords internal
.check_mcmc <- function(y, differencing_option, controls, mcmc_params, regression) {
  .check_data(y, allow_missing = TRUE)
  initial_rows <- if (differencing_option == "D2") 2L else 1L
  # The inherited log-volatility updates require at least three evolution errors.
  if (length(y) < initial_rows + 3L) stop("At least ", initial_rows + 3L, " time points are required for this MCMC sampler.")
  if (sum(!is.na(y)) < 2L || !is.finite(sd(y, na.rm = TRUE)) || sd(y, na.rm = TRUE) <= 0) {
    stop("At least two observed y values with positive sample variance are required for initialisation.")
  }
  for (nm in c("nsave", "nburn", "nskip")) {
    z <- controls[[nm]]
    minimum <- if (nm == "nsave") 1 else 0
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z < minimum || z != floor(z)) {
      stop(nm, " must be an integer >= ", minimum, ".")
    }
  }
  for (nm in setdiff(names(controls), c("nsave", "nburn", "nskip"))) {
    z <- controls[[nm]]
    if (!is.logical(z) || length(z) != 1L || is.na(z)) stop(nm, " must be TRUE or FALSE.")
  }
  if (controls$computeDIC && controls$nsave < 2) stop("computeDIC requires nsave >= 2.")
  total <- controls$nburn + (controls$nskip + 1) * controls$nsave
  if (!is.finite(total) || total > .Machine$integer.max) stop("Requested MCMC iteration count is too large.")
  allowed <- c("mu", "yhat", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean", "rho", if (regression) "beta")
  if (is.list(mcmc_params)) {
    if (!all(vapply(mcmc_params, function(x) is.character(x) && length(x) == 1L, logical(1)))) stop("mcmc_params must contain single character names.")
    mcmc_params <- unlist(mcmc_params, use.names = FALSE)
  }
  if (!is.character(mcmc_params) || anyNA(mcmc_params) || any(!mcmc_params %in% allowed)) {
    stop("mcmc_params must contain only: ", paste(allowed, collapse = ", "), ".")
  }
  unique(mcmc_params)
}

#' Construct the state-evolution operator in time-major order
#'
#' The coefficient vector is c(beta[1, ], beta[2, ], ...). Each row of H
#' produces an initial state or an evolution error. The initial-state priors
#' are proper, so H is square with a unit diagonal and is nonsingular.
#'
#' AR:  beta[t,j] - rho[j] * beta[t-1,j], t >= 2.
#' D1:  beta[t,j] - beta[t-1,j], t >= 2.
#' DAR: beta[2,j] - (1-rho[j]) * beta[1,j] at the boundary;
#'      beta[t,j] - (1+rho[j])*beta[t-1,j] + rho[j]*beta[t-2,j], t >= 3.
#' D2:  beta[t,j] - 2*beta[t-1,j] + beta[t-2,j], t >= 3;
#'      beta[1,j] and beta[2,j] each have an initial-state prior.
#' The DAR boundary follows the original extensions' backfitting convention.
#' It is part of this model specification, not a generic AR prior identity.
#'
#' @param n Number of time points (at least two; at least three for D2).
#' @param p Number of coefficient paths.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar or length-p vector for AR/DAR; ignored for D1/D2.
#' @return Sparse np by np lower-triangular Matrix.
#' @keywords internal
.state_operator <- function(n, p, differencing_option, rho = NULL) {
  differencing_option <- match.arg(differencing_option, c("AR", "D1", "DAR", "D2"))
  minimum <- if (differencing_option == "D2") 3L else 2L
  if (length(n) != 1L || !is.finite(n) || n < minimum || n != floor(n) ||
      length(p) != 1L || !is.finite(p) || p < 1 || p != floor(p)) stop("Invalid state dimensions for ", differencing_option, ".")
  rho <- .rho_vector(rho, p, differencing_option)
  index <- matrix(seq_len(n * p), nrow = n, byrow = TRUE)
  rows <- index[-1, , drop = FALSE]
  lag1 <- switch(differencing_option,
    AR = matrix(rep(-rho, each = n - 1L), nrow = n - 1L),
    D1 = matrix(-1, n - 1L, p),
    DAR = matrix(rep(-(1 + rho), each = n - 1L), nrow = n - 1L),
    D2 = matrix(-2, n - 1L, p))
  if (differencing_option == "DAR") lag1[1, ] <- -(1 - rho)
  if (differencing_option == "D2") lag1[1, ] <- 0
  # Add only the first two subdiagonals in each predictor path: O(np) storage.
  i <- c(seq_len(n * p), as.vector(rows))
  j <- c(seq_len(n * p), as.vector(rows - p))
  x <- c(rep(1, n * p), as.vector(lag1))
  if (differencing_option %in% c("DAR", "D2") && n > 2L) {
    rows2 <- index[-c(1, 2), , drop = FALSE]
    lag2 <- if (differencing_option == "DAR") rep(rho, each = n - 2L) else rep(1, (n - 2L) * p)
    i <- c(i, as.vector(rows2)); j <- c(j, as.vector(rows2 - 2L * p)); x <- c(x, lag2)
  }
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n * p, n * p))
}

#' Compute state-evolution errors
#' @param beta Numeric n by p matrix, or a vector for one path.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar or length-p vector for AR/DAR.
#' @return Matrix of n-1 errors per path, or n-2 for D2.
compute_omega <- function(beta, differencing_option, rho = NULL) {
  beta <- as.matrix(beta)
  n <- nrow(beta); p <- ncol(beta)
  H <- .state_operator(n, p, differencing_option, rho)
  # Multiplying by H enforces the same boundary equations as both samplers.
  states <- matrix(as.numeric(H %*% as.vector(t(beta))), n, p, byrow = TRUE)
  initial_rows <- if (differencing_option == "D2") 2L else 1L
  states[-seq_len(initial_rows), , drop = FALSE]
}

#' Extract the states with initial-state priors
#' @param beta Numeric n by p matrix, or vector for one path.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @return One row for AR/D1/DAR; two rows for D2.
compute_beta0 <- function(beta, differencing_option) {
  differencing_option <- match.arg(differencing_option, c("AR", "D1", "DAR", "D2"))
  beta <- as.matrix(beta)
  initial_rows <- if (differencing_option == "D2") 2L else 1L
  beta[seq_len(initial_rows), , drop = FALSE]
}

#' Draw a scalar truncated normal using log tail probabilities
#'
#' Log probabilities avoid underflow when the untruncated mean lies outside
#' the truncation interval. Moving that mean into the interval would change
#' the distribution, so truncation is applied to the draw only.
#' @param mean Untruncated normal mean.
#' @param sd Strictly positive normal standard deviation.
#' @param lower Finite lower bound.
#' @param upper Finite upper bound, greater than lower.
#' @return One draw in [lower, upper].
#' @keywords internal
.rtrunc_normal <- function(mean, sd, lower, upper) {
  if (any(!is.finite(c(mean, sd, lower, upper))) || sd <= 0 || lower >= upper) stop("Invalid truncated-normal parameters.")
  a <- (lower - mean) / sd; b <- (upper - mean) / sd
  lower_tail <- a <= 0
  logs <- pnorm(c(a, b), lower.tail = lower_tail, log.p = TRUE)
  u <- runif(1)
  terms <- logs + log(c(u, 1 - u))
  largest <- max(terms)
  log_probability <- largest + log(sum(exp(terms - largest)))
  z <- qnorm(log_probability, lower.tail = lower_tail, log.p = TRUE)
  draw <- mean + sd * z
  if (!is.finite(draw)) stop("Truncated-normal draw exceeded floating-point range.")
  # Correct only floating-point roundoff in the final draw, never the mean.
  min(upper, max(lower, draw))
}

#' Sample an AR/DAR autocorrelation parameter from its full conditional
#'
#' Passing all coefficient columns samples one shared rho. Passing one column
#' samples a predictor-specific rho. Both models are weighted regressions in
#' rho, giving a truncated-normal conditional under either supported prior.
#'
#' @param beta Numeric n by p coefficient matrix, or vector for one path.
#'   Set beta and evol_sigma_t2 to NULL to draw from the prior.
#' @param differencing_option Either AR or DAR.
#' @param evol_sigma_t2 Positive finite n by p variance matrix, including initial
#'   states; a length-n vector is accepted for one path. Row one is excluded
#'   from this conditional because its prior is independent of rho.
#' @param prior_rho Increasing length-two vector of bounds within (-1, 1).
#' @param prior_type Either uniform or truncated_normal.
#' @param mu_prior Untruncated mean of the normal prior.
#' @param sigma_prior Positive standard deviation of the normal prior.
#' @param rho_current Retained for call compatibility; the direct conditional
#'   draw no longer needs a previous value or a slice-sampling step.
#' @return One numeric autocorrelation draw within prior_rho.
sample_AR1_param <- function(beta, differencing_option, evol_sigma_t2 = NULL,
                             prior_rho = c(0, 0.99),
                             prior_type = c("uniform", "truncated_normal"),
                             mu_prior = 0.5, sigma_prior = 10,
                             rho_current = NULL) {
  differencing_option <- match.arg(differencing_option, c("AR", "DAR"))
  prior_type <- match.arg(prior_type)
  if (!is.numeric(prior_rho) || length(prior_rho) != 2L || any(!is.finite(prior_rho)) ||
      prior_rho[1] >= prior_rho[2] || any(abs(prior_rho) >= 1)) stop("prior_rho must be increasing bounds strictly inside (-1, 1).")
  if (!is.numeric(mu_prior) || length(mu_prior) != 1L || !is.finite(mu_prior) ||
      !is.numeric(sigma_prior) || length(sigma_prior) != 1L || !is.finite(sigma_prior) || sigma_prior <= 0) stop("mu_prior must be finite and sigma_prior must be positive and finite.")
  prior_draw <- function() {
    if (prior_type == "uniform") runif(1, prior_rho[1], prior_rho[2]) else
      .rtrunc_normal(mu_prior, sigma_prior, prior_rho[1], prior_rho[2])
  }
  if (is.null(beta) && is.null(evol_sigma_t2)) return(prior_draw())
  if (is.null(beta) || is.null(evol_sigma_t2)) stop("Supply both beta and evol_sigma_t2, or set both to NULL.")
  B <- as.matrix(beta); Es <- as.matrix(evol_sigma_t2)
  n <- nrow(B)
  if (!is.numeric(B) || any(!is.finite(B)) || n < 2L || ncol(B) < 1L ||
      !identical(dim(B), dim(Es))) stop("beta and evol_sigma_t2 must have matching n by p dimensions, n >= 2.")
  .check_variances(Es, "evol_sigma_t2")
  if (differencing_option == "AR") {
    response <- B[-1, , drop = FALSE]
    predictor <- B[-n, , drop = FALSE]
  } else {
    # The first DAR residual is beta[2] - beta[1] + rho*beta[1].
    response <- diff(B)
    predictor <- rbind(-B[1, , drop = FALSE], if (n > 2L) diff(B)[- (n - 1L), , drop = FALSE])
  }
  # Standardising by SDs avoids an explicit matrix inverse or clipped variances.
  x <- predictor / sqrt(Es[-1, , drop = FALSE])
  y <- response / sqrt(Es[-1, , drop = FALSE])
  precision <- sum(x^2)
  linear <- sum(x * y)
  if (prior_type == "truncated_normal") {
    precision <- precision + 1 / sigma_prior^2
    linear <- linear + mu_prior / sigma_prior^2
  }
  if (!is.finite(precision) || !is.finite(linear)) stop("rho conditional exceeded floating-point range; check scales.")
  if (precision == 0) return(prior_draw())
  .rtrunc_normal(linear / precision, 1 / sqrt(precision), prior_rho[1], prior_rho[2])
}

#' Build the evolution-prior precision from the state equations
#' @param evol_sigma_t2 Positive finite n by p matrix, including initial states.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar or length-p vector for AR/DAR.
#' @return Sparse symmetric np by np Matrix H' diag(1 / variance) H.
#' @keywords internal
.evolution_precision <- function(evol_sigma_t2, differencing_option, rho = NULL) {
  evol_sigma_t2 <- as.matrix(evol_sigma_t2)
  .check_variances(evol_sigma_t2, "evol_sigma_t2")
  H <- .state_operator(nrow(evol_sigma_t2), ncol(evol_sigma_t2), differencing_option, rho)
  weighted_H <- Matrix::Diagonal(x = 1 / sqrt(as.vector(t(evol_sigma_t2)))) %*% H
  Matrix::crossprod(weighted_H)
}

#' Build the joint regression posterior precision
#' @param obs_sigma_t2 Positive finite length-n observation variance vector.
#' @param evol_sigma_t2 Positive finite n by p evolution variance matrix.
#' @param XtX Block diagonal np by np Matrix from build_XtX(X).
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar or length-p vector for AR/DAR.
#' @return Sparse symmetric np by np posterior precision.
#' @keywords internal
.reg_precision <- function(obs_sigma_t2, evol_sigma_t2, XtX, differencing_option, rho = NULL) {
  evol_sigma_t2 <- as.matrix(evol_sigma_t2)
  n <- nrow(evol_sigma_t2); p <- ncol(evol_sigma_t2)
  .check_variances(obs_sigma_t2, "obs_sigma_t2")
  if (length(obs_sigma_t2) != n || !identical(as.integer(dim(XtX)), as.integer(c(n*p, n*p)))) stop("Variance and XtX dimensions do not match.")
  Qevol <- .evolution_precision(evol_sigma_t2, differencing_option, rho)
  # XtX contains x[t,] x[t,]' blocks, not the usual p by p crossprod(X).
  Qobs <- Matrix::Diagonal(x = rep(1 / as.numeric(obs_sigma_t2), each = p)) %*% XtX
  Matrix::forceSymmetric(Qobs + Qevol)
}

#' Equilibrate a Gaussian precision without changing its distribution
#' @param Q Sparse symmetric positive-definite precision matrix.
#' @return List containing the scaled precision and diagonal scale vector.
#' @keywords internal
.scale_precision <- function(Q) {
  diagonal <- Matrix::diag(Q)
  if (any(!is.finite(Q@x)) || any(!is.finite(diagonal)) || any(diagonal <= 0)) stop("Posterior precision contains non-finite values or a non-positive diagonal; check data and variance scales.")
  scale <- 1 / sqrt(diagonal)
  S <- Matrix::Diagonal(x = scale)
  list(Q = Matrix::forceSymmetric(S %*% Q %*% S), scale = scale)
}

#' Convert a sparse precision to spam without deprecated symmetric coercions
#' @param Q Sparse symmetric Matrix.
#' @return A spam matrix with both triangles represented.
#' @keywords internal
.as_spam <- function(Q) {
  spam::as.spam.dgCMatrix(methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix"))
}

#' Sample a Gaussian using its canonical parameters
#' @param Q Sparse positive-definite precision matrix.
#' @param b Finite linear term in exp(-x'Qx/2 + b'x).
#' @param chol0 Optional spam symbolic Cholesky structure for this ordering.
#' @return Numeric draw with mean solve(Q, b) and covariance solve(Q).
#' @keywords internal
.sample_gaussian_precision <- function(Q, b, chol0 = NULL) {
  scaled <- .scale_precision(Q)
  b <- as.numeric(b)
  if (length(b) != nrow(Q) || any(!is.finite(b))) stop("Gaussian linear term has invalid dimensions or non-finite values.")
  b_scaled <- scaled$scale * b
  if (!is.null(chol0)) {
    # A cached symbolic pattern can become obsolete when rho moves from zero.
    # Retry with Matrix if spam cannot update it; neither path adds a ridge.
    draw <- tryCatch(as.numeric(spam::rmvnorm.canonical(
      n = 1, b = b_scaled, Q = .as_spam(scaled$Q), Rstruct = chol0)),
      error = function(e) {
        warning("Cached spam factorisation failed; retrying Matrix: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    if (!is.null(draw)) return(scaled$scale * draw)
  }
  # For Qs = R'R, R^{-1}(R'^{-1}b + z) has the required mean and covariance.
  # Explicitly disable permutation because these triangular solves use b's order.
  R <- tryCatch(Matrix::chol(scaled$Q, pivot = FALSE), error = function(e) {
    stop("Cholesky failed after diagonal scaling. Check variance ranges and predictor scales; regression may benefit from use_backfitting = TRUE. Original error: ", conditionMessage(e), call. = FALSE)
  })
  draw <- Matrix::solve(R, Matrix::solve(Matrix::t(R), b_scaled) + rnorm(length(b)))
  scaled$scale * as.numeric(draw)
}

#' Initialise a sparse Cholesky structure for joint TVP regression
#' @param obs_sigma_t2 Positive finite length-n observation variances.
#' @param evol_sigma_t2 Positive finite n by p evolution variances.
#' @param XtX Block diagonal np by np Matrix from build_XtX(X).
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar or length-p vector for AR/DAR.
#' @return A spam.chol.NgPeyton factor for reuse with sampleBTF_reg().
#' @note Reuse avoids symbolic analysis when the sparsity pattern is unchanged;
#' numeric values are factorised again at each draw.
initCholReg.spam <- function(obs_sigma_t2, evol_sigma_t2, XtX,
                             differencing_option = "D1", rho = NULL) {
  Q <- .reg_precision(obs_sigma_t2, evol_sigma_t2, XtX, differencing_option, rho)
  spam::chol.spam(.as_spam(.scale_precision(Q)$Q))
}

#' Draw all time-varying regression coefficients jointly
#'
#' The observation model is y[t] = sum(X[t, ] * beta[t, ]) + error[t].
#' Combining its Gaussian likelihood with H beta gives precision
#' Q = Z' diag(1 / obs_sigma_t2) Z + H' diag(1 / evol_sigma_t2) H,
#' where row t of Z contains X[t, ] in the t-th time block. Sparse Cholesky
#' and triangular solves produce a draw without forming the covariance inverse.
#'
#' @param y Finite numeric length-n response vector; no missing values.
#' @param X Finite numeric n by p predictor matrix; include an intercept explicitly.
#' @param obs_sigma_t2 Positive finite length-n observation variances.
#' @param evol_sigma_t2 Positive finite n by p evolution variances, including
#'   the initial-state variances in row one (rows one and two for D2).
#' @param XtX Block diagonal np by np Matrix from build_XtX(X) for the same X.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar shared rho or length-p predictor-specific vector for AR/DAR.
#' @param chol0 Optional spam Cholesky structure from initCholReg.spam(); NULL
#'   selects Matrix. A failed cached update retries with Matrix.
#' @return Numeric n by p matrix of sampled coefficient paths.
#' @note See .state_operator() for the boundary equations. btf_reg() handles
#' missing responses before calling this conditional sampler.
sampleBTF_reg <- function(y, X, obs_sigma_t2, evol_sigma_t2, XtX,
                          differencing_option = "D1", rho = NULL, chol0 = NULL) {
  .check_data(y, X)
  n <- nrow(X); p <- ncol(X)
  if (!identical(dim(as.matrix(evol_sigma_t2)), dim(X))) stop("evol_sigma_t2 must have the same dimensions as X.")
  Q <- .reg_precision(obs_sigma_t2, evol_sigma_t2, XtX, differencing_option, rho)
  # Flatten by time, matching H and the block-diagonal likelihood precision.
  b <- as.vector(t(X * as.numeric(y / obs_sigma_t2)))
  matrix(.sample_gaussian_precision(Q, b, chol0), n, p, byrow = TRUE)
}

#' Update time-varying coefficients by a random-order backfitting sweep
#'
#' Each predictor path is drawn conditional on the current other paths. A
#' sweep leaves the joint conditional distribution invariant, but successive
#' sweeps can be correlated. Each path uses a banded n by n Cholesky solve.
#' With a running residual, work per sweep is O(np) for these state operators.
#'
#' @param y Finite numeric length-n response vector; no missing values.
#' @param X Finite numeric n by p predictor matrix.
#' @param beta Finite numeric n by p matrix of current coefficient paths.
#' @param obs_sigma_t2 Positive finite length-n observation variances.
#' @param evol_sigma_t2 Positive finite n by p evolution variances.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar shared rho or length-p vector for AR/DAR.
#' @return Numeric n by p coefficient matrix after one complete sweep.
#' @note Zero predictor values contribute zero likelihood precision; no
#' division of the response by a predictor is necessary.
sampleBTF_reg_backfit <- function(y, X, beta, obs_sigma_t2, evol_sigma_t2,
                                  differencing_option = "D1", rho = NULL) {
  .check_data(y, X)
  n <- nrow(X); p <- ncol(X)
  if (!is.numeric(beta) || !identical(dim(beta), dim(X)) || any(!is.finite(beta)) ||
      !identical(dim(as.matrix(evol_sigma_t2)), dim(X))) stop("beta and evol_sigma_t2 must match X's dimensions; beta must be finite.")
  .check_variances(obs_sigma_t2, "obs_sigma_t2")
  if (length(obs_sigma_t2) != n) stop("obs_sigma_t2 must have length n.")
  obs_sigma_t2 <- as.numeric(obs_sigma_t2)
  differencing_option <- match.arg(differencing_option, c("AR", "D1", "DAR", "D2"))
  rho <- .rho_vector(rho, p, differencing_option)
  residual <- y - rowSums(X * beta)
  for (j in sample.int(p)) {
    # Add back the old j-th contribution, sample it, then remove the new one.
    partial_residual <- residual + X[, j] * beta[, j]
    Q <- .evolution_precision(evol_sigma_t2[, j, drop = FALSE],
                              differencing_option, if (is.null(rho)) NULL else rho[j]) +
      Matrix::Diagonal(x = X[, j]^2 / obs_sigma_t2)
    b <- partial_residual * X[, j] / obs_sigma_t2
    beta[, j] <- .sample_gaussian_precision(Q, b)
    residual <- partial_residual - X[, j] * beta[, j]
  }
  beta
}

#' Build the univariate Gaussian posterior precision
#' @param obs_sigma_t2 Positive length-n observation variances; Inf gives zero
#'   likelihood precision, useful for an absent observation contribution.
#' @param evol_sigma_t2 Positive finite length-n evolution variances, including
#'   one initial-state variance (two for D2).
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar autocorrelation for AR/DAR.
#' @return Sparse symmetric banded n by n Matrix.
build_Q <- function(obs_sigma_t2, evol_sigma_t2, differencing_option = "D1", rho = NULL) {
  .check_variances(obs_sigma_t2, "obs_sigma_t2", allow_inf = TRUE)
  if (length(obs_sigma_t2) != length(evol_sigma_t2)) stop("Observation and evolution variances must have the same length.")
  .evolution_precision(matrix(evol_sigma_t2, ncol = 1), differencing_option, rho) +
    Matrix::Diagonal(x = 1 / as.numeric(obs_sigma_t2))
}

# ============================================================================
# Shrinkage priors, log-volatility updates, and MCMC utilities
# ============================================================================

#----------------------------------------------------------------------------
#' Compute X'X
#'
#' Build the \code{np x np} matrix XtX using the Matrix package
#' @param X \code{n x p} matrix of predictors
#' @return Block diagonal \code{np x np} Matrix (object) where each \code{p x p} block is \code{tcrossprod(matrix(X[t,]))}
#'
#' @note X'X is a one-time computing cost. Special cases may have more efficient computing options,
#' but the Matrix representation is important for efficient computations within the sampler.
#'
#' @import Matrix
#' @export
build_XtX <- function(X){

  # Store the dimensions:
  n = nrow(X); p = ncol(X)

  # Store the matrix
  XtX = Matrix::bandSparse(n*p, k = 0, diagonals = list(rep(1,n*p)), symmetric = TRUE)

  t.seq.p = seq(1, n*(p+1), by = p)

  for(t in 1:n){
    t.ind = t.seq.p[t]:(t.seq.p[t+1]-1)
    XtX[t.ind, t.ind] = tcrossprod(matrix(X[t,]))
  }
  XtX
}

#----------------------------------------------------------------------------
#' Initialise the evolution error variance parameters
#'
#' Compute initial values for evolution error variance parameters under the various options:
#' dynamic horseshoe prior ('DHS'), horseshoe prior ('HS'),
#' Bayesian lasso ('BL'), normal stochastic volatility ('SV'),
#' or normal-inverse-gamma prior ('NIG').
#'
#' @param omega \code{n x p} matrix of evolution errors
#' @param evol_error the evolution error distribution; must be one of
#' DHS (dynamic horseshoe), HS (horseshoe), BL (Bayesian lasso), SV (stochastic volatility), or NIG (normal-inverse-gamma)
#' @return List of relevant components: \code{sigma_wt}, the \code{n x p} matrix of evolution standard deviations,
#' and additional parameters associated with the DHS and HS priors.
#' @export
initEvolParams <- function(omega, evol_error = "DHS"){

  # Check:
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") ||(evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  if(evol_error == "DHS") return(initDHS(omega))

  if(evol_error == "HS"){
    tauLambdaj = 1/omega^2;
    xiLambdaj = 1/(2*tauLambdaj); tauLambda = 1/(2*colMeans(xiLambdaj)); xiLambda = 1/(tauLambda + 1)

    # Parameters to store/return:
    return(list(sigma_wt = 1/sqrt(tauLambdaj), tauLambdaj = tauLambdaj, xiLambdaj = xiLambdaj, tauLambda = tauLambda, xiLambda = xiLambda))
  }
  if(evol_error == "BL"){
    tau_j = abs(omega); lambda2 = mean(tau_j)
    return(list(sigma_wt = tau_j, tau_j = tau_j, lambda2 = lambda2))
  }
  if(evol_error == "SV") return(initSV(omega))
  if(evol_error == "NIG") return(list(sigma_wt = tcrossprod(rep(1,n), apply(omega, 2, function(x) sd(x, na.rm=TRUE)))))
}

#----------------------------------------------------------------------------
#' Initialise the evolution error variance parameters
#'
#' Compute initial values for evolution error variance parameters under the dynamic horseshoe prior
#'
#' @param omega \code{n x p} matrix of evolution errors
#' @return List of relevant components: the \code{n x p} evolution error SD \code{sigma_wt},
#' the \code{n x p} log-volatility \code{ht}, the \code{p x 1} log-vol unconditional mean(s) \code{dhs_mean},
#' the \code{p x 1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' the \code{n x p} log-vol innovation SD \code{sigma_eta_t} from the PG priors,
#' the \code{p x 1} initial log-vol SD \code{sigma_eta_0},
#' and the mean of log-vol means \code{dhs_mean0} (relevant when \code{p > 1})
#' @export
initDHS <- function(omega){

  # "Local" number of time points
  omega = as.matrix(omega)
  n = nrow(omega); p = ncol(omega)

  # Initialise the log-volatilities:
  ht = log(omega^2 + 0.0001)

  # Initialise the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = apply(ht, 2, function(x){
    params = try(arima(x, c(1,0,0))$coef, silent = TRUE); if(paste(class(params)) == "try-error") params = c(0.8, mean(x)/(1 - 0.8))
    params
  })
  dhs_mean = arCoefs[2,]; dhs_phi = arCoefs[1,]; dhs_mean0 = mean(dhs_mean)

  # Initialise the SD of log-vol innovations simply using the expectation:
  sigma_eta_t = matrix(pi, nrow = n-1, ncol = p)
  sigma_eta_0 = rep(pi, p) # Initial value

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}

#----------------------------------------------------------------------------
#' Initialise the parameters for the initial state variance
#'
#' The initial state SDs are assumed to follow half-Cauchy priors, C+(0,A),
#' where the SDs may be common or distinct among the states.
#'
#' This function initalises the parameters for a PX-Gibbs sampler.
#'
#' @param mu0 Vector or matrix of initial states (one row per initial time);
#'   each entry receives a scale when commonSD is FALSE.
#' @param commonSD logical; if TRUE, use common SDs (otherwise distinct)
#' @return List of relevant components:
#' the \code{p x 1} evolution error SD \code{sigma_w0},
#' the \code{p x 1} parameter-expanded RV's \code{px_sigma_w0},
#' and the corresponding global scale parameters
#' \code{sigma_00} and \code{px_sigma_00} (ignore if commonSD)
#' @export
initEvol0 <- function(mu0, commonSD = TRUE){

  p = length(mu0)

  # Common or distinct:
  if(commonSD) {
    sigma_w0 = rep(mean(abs(mu0)), p)
  } else  sigma_w0 = abs(mu0)

  # Initialise at 1 for simplicity:
  px_sigma_w0 = rep(1, p)

  sigma_00 = px_sigma_00 = 1

  list(sigma_w0 = sigma_w0, px_sigma_w0 = px_sigma_w0, sigma_00 = sigma_00, px_sigma_00 = px_sigma_00)
}

#----------------------------------------------------------------------------
#' Initialise the stochastic volatility parameters
#'
#' Compute initial values for normal stochastic volatility parameters.
#' The model assumes an AR(1) for the log-volatility.
#'
#' @param omega \code{n x p} matrix of errors
#' @return List of relevant components: \code{sigma_wt}, the \code{n x p} matrix of standard deviations,
#' and additional parameters (unconditional mean, AR(1) coefficient, and standard deviation).
#' @export
initSV <- function(omega){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  # log-volatility:
  ht = log(omega^2 + 0.0001)

  # AR(1) pararmeters: check for error in initialization too
  svParams = apply(ht, 2, function(x){
    ar_fit = try(arima(x, c(1,0,0)), silent = TRUE)
    if(paste(class(ar_fit)) != "try-error") {
      params = c(ar_fit$coef[2], ar_fit$coef[1], sqrt(ar_fit$sigma2))
    } else params = c(mean(x)/(1 - 0.8),0.8, 1)
    params
  }); rownames(svParams) = c("intercept", "ar1", "sig")

  # SDs, log-vols, and other parameters:
  return(list(sigma_wt = exp(ht/2), ht = ht, svParams = svParams))
}

#----------------------------------------------------------------------------
#' Sample the parameters for the initial state variance
#'
#' The initial state SDs are assumed to follow half-Cauchy priors, C+(0,A),
#' where the SDs may be common or distinct among the states.
#'
#' This function samples the parameters for a PX-Gibbs sampler.
#'
#' @param mu0 Vector or matrix of initial states (one row per initial time);
#'   each entry receives a scale when commonSD is FALSE.
#' @param evolParams0 list of relevant components (see below)
#' @param commonSD logical; if TRUE, use common SDs (otherwise distinct)
#' @param A prior scale parameter from the half-Cauchy prior, C+(0,A)
#' @return List of relevant components:
#' the \code{p x 1} evolution error SD \code{sigma_w0}
#' and the \code{p x 1} parameter-expanded RV's \code{px_sigma_w0}
#' @export
sampleEvol0 <- function(mu0, evolParams0, commonSD = FALSE, A = 1){

  # Store length locally:
  p = length(mu0)

  # For numerical stability:
  mu02offset = any(mu0^2 < 10^-16)*max(10^-8, mad(mu0)/10^6)
  mu02 = mu0^2 + mu02offset

  if(commonSD){
    # (Common) standard deviations:
    evolParams0$sigma_w0 = rep(1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(mu02)/2 + evolParams0$px_sigma_w0[1])), p)

    # (Common) paramater expansion:
    evolParams0$px_sigma_w0 = rep(rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0[1]^2 + 1/A^2), p)

  } else {
    # (Distinct) standard deviations:
    evolParams0$sigma_w0 = 1/sqrt(rgamma(n = p, shape = 1/2 + 1/2, rate = mu02/2 + evolParams0$px_sigma_w0))

    # (distinct) paramater expansion:
    evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/evolParams0$sigma_00^2)

    # Global standard deviations:
    evolParams0$sigma_00 = 1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(evolParams0$px_sigma_w0) + evolParams0$px_sigma_00))

    # (Global) parameter expansion:
    evolParams0$px_sigma_00 = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_00^2 + 1/A^2)
  }

  # And return the list:
  evolParams0
}

#----------------------------------------------------------------------------
#' Sample evolution error variance parameters
#'
#' Compute one draw of evolution error variance parameters under the various options:
#' \itemize{
#' \item dynamic horseshoe prior ('DHS');
#' \item horseshoe prior ('HS');
#' \item Bayesian lasso ('BL');
#' \item normal stochastic volatility ('SV');
#' \item normal-inverse-gamma prior ('NIG').
#' }
#'
#' @param omega \code{n x p} matrix of evolution errors
#' @param evolParams list of parameters pertaining to each \code{evol_error} type to be updated
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @param evol_error the evolution error distribution; must be one of
#' DHS (dynamic horseshoe), HS (horseshoe), BL (Bayesian lasso), SV (stochastic volatility), or NIG (normal-inverse-gamma)
#' @return List of relevant components in \code{evolParams}: \code{sigma_wt}, the \code{n x p} matrix of evolution standard deviations,
#' and additional parameters associated with the DHS and HS priors.
#'
#' @note The list \code{evolParams} is specific to each \code{evol_error} type,
#' but in each case contains the evolution error standard deviations \code{sigma_wt}.
#'
#' @note sigma_e scales the DHS and BL priors in this implementation.
#' The retained HS, SV and NIG updates do not use it. Use sigma_e = 1 to
#' omit observation-scale adjustment in the branches that support it.
#'
#' @import stochvol
#' @export
sampleEvolParams <- function(omega, evolParams,  sigma_e = 1, evol_error = "DHS"){

  # Check:
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") || (evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  if(evol_error == "DHS") return(sampleDSP(omega, evolParams, sigma_e))

  if(evol_error == "HS"){

    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset

    # Local scale params:
    evolParams$tauLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$xiLambdaj + hsInput2/2), nrow = n)
    evolParams$xiLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$tauLambdaj + tcrossprod(rep(1,n), evolParams$tauLambda)), nrow = n)

    # Global scale params:
    evolParams$tauLambda = rgamma(n = p, shape = 0.5 + n/2, colSums(evolParams$xiLambdaj) + evolParams$xiLambda)

    evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1)

    evolParams$sigma_wt = 1/sqrt(evolParams$tauLambdaj)

    return(evolParams)
  }
  if(evol_error == "BL"){

    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset

    # 1/tau_j^2 is inverse-gaussian (NOTE: this is very slow!)
    evolParams$tau_j = matrix(sapply(matrix(hsInput2), function(x){1/sqrt(rig(n = 1,
                                            mean = sqrt(evolParams$lambda2*sigma_e^2/x), # already square the input
                                            scale = 1/evolParams$lambda2))}), nrow = n)
    # Note: should be better priors for lambda2
    evolParams$lambda2 = rgamma(n = 1,
                                shape = 1 + n*p,
                                rate = 2 + sum(evolParams$tau_j^2)/2)

    # For Bayesian lasso, scale by sigma_e:
    evolParams$sigma_wt = sigma_e*evolParams$tau_j

    return(evolParams)
  }
  if(evol_error == "SV") return(sampleSVparams(omega = omega, svParams = evolParams))
  if(evol_error == "NIG") {
    evolParams = list(sigma_wt = tcrossprod(rep(1,n),
                                            apply(omega, 2,
                                                  function(x) 1/sqrt(rgamma(n = 1, shape = n/2 + 0.01, rate = sum(x^2)/2 + 0.01)))))
    return(evolParams)
  }
}

#----------------------------------------------------------------------------
#' Sample the dynamic shrinkage process parameters
#'
#' Compute one draw for each of the parameters in the dynamic shrinkage process
#' for the special case in which the shrinkage parameter \code{kappa ~ Beta(alpha, beta)}
#' with \code{alpha = beta}. The primary example is the dynamic horseshoe process with
#' \code{alpha = beta = 1/2}.
#'
#' @param omega \code{n x p} matrix of evolution errors
#' @param evolParams list of parameters to be updated (see Value below)
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @param prior_dhs_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#' @param alphaPlusBeta For the symmetric prior kappa ~ Beta(alpha, beta) with alpha=beta,
#' specify the sum [alpha + beta]
#' @return List of relevant components:
#' \itemize{
#' \item the \code{n x p} evolution error standard deviations \code{sigma_wt},
#' \item the \code{n x p} log-volatility \code{ht}, the \code{p x 1} log-vol unconditional mean(s) \code{dhs_mean},
#' \item the \code{p x 1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' \item the \code{n x p} log-vol innovation standard deviations \code{sigma_eta_t} from the Polya-Gamma priors,
#' \item the \code{p x 1} initial log-vol SD \code{sigma_eta_0},
#' \item and the mean of log-vol means \code{dhs_mean0} (relevant when \code{p > 1})
#' }
#'
#' @note The priors induced by \code{prior_dhs_phi} all imply a stationary (log-) volatility process.
#'
#' @import BayesLogit
#' @export
sampleDSP <- function(omega, evolParams, sigma_e = 1, prior_dhs_phi = c(10,2), alphaPlusBeta = 1){

  # Store the DSP parameters locally:
  ht = evolParams$ht; dhs_mean = evolParams$dhs_mean; dhs_phi = evolParams$dhs_phi; sigma_eta_t = evolParams$sigma_eta_t; sigma_eta_0 = evolParams$sigma_eta_0; dhs_mean0 = evolParams$dhs_mean0

  # "Local" number of time points
  ht = as.matrix(ht)
  n = nrow(ht); p = ncol(ht)

  # Sample the log-volatilities using AWOL sampler
  ht = sampleLogVols(h_y = omega, h_prev = ht, h_mu = dhs_mean, h_phi=dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0)

  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - tcrossprod(rep(1,n), dhs_mean)

  # Sample AR(1) parameters
    # Note: dhs_phi = 0 means non-dynamic HS, while dhs_phi = 1 means RW, in which case we don't sample either
  if(!all(dhs_phi == 0) && !all(dhs_phi == 1)) dhs_phi = sampleAR1(h_yc = ht_tilde, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, prior_dhs_phi = prior_dhs_phi)

  # Sample the evolution error SD of log-vol (i.e., Polya-Gamma mixing weights)
  eta_t = ht_tilde[-1,] - tcrossprod(rep(1,n-1), dhs_phi)*ht_tilde[-n, ]       # Residuals
  sigma_eta_t = matrix(1/sqrt(rpg(num = (n-1)*p, h = alphaPlusBeta, z = eta_t)), ncol = p) # Sample
  sigma_eta_0 = 1/sqrt(rpg(num = p, h = 1, z = ht_tilde[1,]))                # Sample the inital

  # Sample the unconditional mean(s), unless dhs_phi = 1 (not defined)
  if(!all(dhs_phi == 1)){
    if(p > 1){
      # Assume a hierarchy of the global shrinkage params across j=1,...,p
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = dhs_mean0);
      dhs_mean = muSample$dhs_mean
      dhs_mean0 = sampleLogVolMu0(h_mu = dhs_mean, h_mu0 = dhs_mean0, dhs_mean_prec_j = muSample$dhs_mean_prec_j, h_log_scale = log(sigma_e^2))
    } else {
      # p = 1
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = log(sigma_e^2));
      dhs_mean = dhs_mean0 = muSample$dhs_mean # save dhs_mean0 = dhs_mean for coding convenience later
    }
  } else {dhs_mean = rep(0, p); dhs_mean0 = 0} # When RW for log-vols, fix unconditional mean for identifiability

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  # Return the same list, but with the new values
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}

#----------------------------------------------------------------------------
#' Sample the AR(1) unconditional means
#'
#' Compute one draw of the unconditional means in an AR(1) model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h the \code{n x p} matrix of log-volatilities
#' @param h_mu the \code{p x 1} vector of previous means
#' @param h_phi the \code{p x 1} vector of AR(1) coefficient(s)
#' @param h_sigma_eta_t the \code{n x p} matrix of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the standard deviations of initial log-vols
#' @param h_log_scale prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2) or dhs_mean0
#'
#' @return a list containing
#' \itemize{
#' \item the sampled mean(s) \code{dhs_mean} and
#' \item the sampled precision(s) \code{dhs_mean_prec_j} from the Polya-Gamma parameter expansion
#'}
#'
#' @import BayesLogit
#' @export
sampleLogVolMu <- function(h, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0, h_log_scale = 0){

  # Compute "local" dimensions:
  n = nrow(h); p = ncol(h)

  # Sample the precision term(s)
  dhs_mean_prec_j = rpg(num = p, h = 1, z = h_mu - h_log_scale)

  # Now, form the "y" and "x" terms in the (auto)regression
  y_mu = (h[-1,] - tcrossprod(rep(1,n-1), h_phi)*h[-n,])/h_sigma_eta_t;
  x_mu = tcrossprod(rep(1,n-1), 1 - h_phi)/h_sigma_eta_t

  # Include the initial sd
  y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu);
  x_mu = rbind(1/h_sigma_eta_0, x_mu)

  # Posterior SD and mean:
  postSD = 1/sqrt(colSums(x_mu^2) + dhs_mean_prec_j)
  postMean = (colSums(x_mu*y_mu) + h_log_scale*dhs_mean_prec_j)*postSD^2
  dhs_mean = rnorm(n = p, mean = postMean, sd = postSD)

  list(dhs_mean = dhs_mean, dhs_mean_prec_j = dhs_mean_prec_j)
}

#----------------------------------------------------------------------------
#' Sample the mean of AR(1) unconditional means
#'
#' Compute one draw of the mean of unconditional means in an AR(1) model with Gaussian innovations
#' and time-dependent innovation variances (for p > 1). More generally, the sampler
#' applies to the "mean" parameter (on the log-scale) for a Polya-Gamma parameter expanded
#' hierarchical model.
#'
#' @param h_mu the \code{p x 1} vector of means
#' @param h_mu0 the previous mean of unconditional means
#' @param dhs_mean_prec_j the \code{p x 1} vector of precisions (from the Polya-Gamma parameter expansion)
#' @param h_log_scale prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2)
#'
#' @return The sampled mean parameter \code{dhs_mean0}
#'
#' @note This sampler is particularly for \code{p > 1} and the setting in which we want hierarchical
#' shrinkage effects, e.g. predictor- and time-dependent shrinkage, predictor-dependent shrinkage,
#' and global shrinkage, with a natural hierarchical ordering.
#'
#' @import BayesLogit
#' @export
sampleLogVolMu0 <- function(h_mu, h_mu0, dhs_mean_prec_j, h_log_scale = 0){

  dhs_mean_prec_0 = rpg(num = 1, h = 1, z = h_mu0 - h_log_scale)

  # Sample the common mean parameter:
  postSD = 1/sqrt(sum(dhs_mean_prec_j) + dhs_mean_prec_0)
  postMean = (sum(dhs_mean_prec_j*h_mu) + dhs_mean_prec_0*h_log_scale)*postSD^2
  rnorm(n = 1, mean = postMean, sd = postSD)
}

#----------------------------------------------------------------------------
#' Sample the AR(1) coefficient(s)
#'
#' Compute one draw of the AR(1) coefficient in a model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h_yc the \code{n x p} matrix of centered log-volatilities
#' (i.e., the log-vols minus the unconditional means \code{dhs_mean})
#' @param h_phi the \code{p x 1} vector of previous AR(1) coefficient(s)
#' @param h_sigma_eta_t the \code{n x p} matrix of log-vol innovation standard deviations
#' @param prior_dhs_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#'
#' @return \code{p x 1} vector of sampled AR(1) coefficient(s)
#'
#' @note For the standard AR(1) case, \code{p = 1}. However, the function applies more
#' generally for sampling \code{p > 1} independent AR(1) processes (jointly).
#'
#' @import truncdist
#' @export
sampleAR1 <- function(h_yc, h_phi, h_sigma_eta_t, prior_dhs_phi = NULL){

  # Compute dimensions:
  n = nrow(h_yc); p = ncol(h_yc)

  # Loop over the j=1:p
  for(j in 1:p){

    # Compute "regression" terms for dhs_phi_j:
    y_ar = h_yc[-1,j]/h_sigma_eta_t[,j] # Standardized "response"
    x_ar = h_yc[-n,j]/h_sigma_eta_t[,j] # Standardized "predictor"

    # Using Beta distribution:
    if(!is.null(prior_dhs_phi)){

      # Check to make sure the prior params make sense
      if(length(prior_dhs_phi) != 2) stop('prior_dhs_phi must be a numeric vector of length 2')

      dhs_phi01 = (h_phi[j] + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])

      # Slice sampler when using Beta prior:
      dhs_phi01 = uni.slice(dhs_phi01, g = function(x){
        -0.5*sum((y_ar - (2*x - 1)*x_ar)^2) +
          dbeta(x, shape1 = prior_dhs_phi[1], shape2 = prior_dhs_phi[2], log = TRUE)
      }, lower = 0, upper = 1)[1]#}, lower = 0.005, upper = 0.995)[1] #

      h_phi[j] = 2*dhs_phi01 - 1

    } else {
      # For h_phi ~ Unif(-1, 1), the posterior is truncated normal
      h_phi[j] = rtrunc(n = 1, spec = 'norm',
                        a = -1, b = 1,
                        mean = sum(y_ar*x_ar)/sum(x_ar^2),
                        sd = 1/sqrt(sum(x_ar^2)))
    }
  }
  h_phi
}

#----------------------------------------------------------------------------
#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.

  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }

  # Find the log density at the initial point, if not already known.

  if (is.null(gx0))
  {
  gx0 <- g(x0)
  }

  # Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      if (g(L)<=logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      if (g(R)<=logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }

  # Shrink interval to lower and upper bounds.
  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.
  repeat
  {
    x1 <- runif(1,L,R)

    gx1 <- g(x1)

    if (gx1>=logy) break

    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.
  attr(x1,"log.density") <- gx1
  return (x1)

}

#----------------------------------------------------------------------------
#' Sample from an inverse-Gaussian distribution
#'
#' Using code from the \code{mgcv} package
#'
#' @param n the number of deviates required. If this has length > 1 then the length is taken as the number of deviates required.
#' @param mean vector of mean values.
#' @param scale vector of scale parameter values (lambda)
#' @return Numeric vector of inverse-Gaussian draws.
rig <- function (n, mean, scale)
{
  if (length(n) > 1)
    n <- length(n)
  x <- y <- rnorm(n)^2
  mys <- mean * scale * y
  mu <- 0 * y + mean
  mu2 <- mu^2
  ind <- mys < .Machine$double.eps^-0.5
  x[ind] <- mu[ind] * (1 + 0.5 * (mys[ind] - sqrt(mys[ind] *
                                                    4 + mys[ind]^2)))
  x[!ind] <- mu[!ind]/mys[!ind]
  ind <- runif(n) > mean/(mean + x)
  x[ind] <- mu2[ind]/x[ind]
  x
}

#----------------------------------------------------------------------------
#' Sample the latent log-volatilities
#'
#' Compute one draw of the log-volatilities using a discrete mixture of Gaussians
#' approximation to the likelihood (see Omori, Chib, Shephard, and Nakajima, 2007)
#' where the log-vols are assumed to follow an AR(1) model with time-dependent
#' innovation variances. More generally, the code operates for \code{p} independent
#' AR(1) log-vol processes to produce an efficient joint sampler in \code{O(np)} time.
#'
#' @param h_y the \code{n x p} matrix of data, which follow independent SV models
#' @param h_prev the \code{n x p} matrix of the previous log-vols
#' @param h_mu the \code{p x 1} vector of log-vol unconditional means
#' @param h_phi the \code{p x 1} vector of log-vol AR(1) coefficients
#' @param h_sigma_eta_t the \code{n x p} matrix of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the \code{p x 1} vector of initial log-vol innovation standard deviations
#'
#' @return \code{n x p} matrix of simulated log-vols
#'
#' @note For Bayesian trend filtering, \code{p = 1}. More generally, the sampler allows for
#' \code{p > 1} but assumes (contemporaneous) independence across the log-vols for \code{j = 1,...,p}.
#'
#' @import Matrix
#' @import BayesLogit
sampleLogVols <- function(h_y, h_prev, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0){

  # Compute dimensions:
  h_prev = as.matrix(h_prev) # Just to be sure (n x p)
  n = nrow(h_prev); p = ncol(h_prev)

  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)

  # Add an offset: common for all times, but distinct for each j=1,...,p
  yoffset = tcrossprod(rep(1,n),
                       apply(as.matrix(h_y), 2,
                             function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))

  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)

  # Sample the mixture components
  z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)

  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = matrix(m_st[z], nrow=n); v_st2_all = matrix(v_st2[z], nrow=n)

  # Joint AWOL sampler for j=1,...,p:

  # Constant (but j-specific) mean
  h_mu_all = tcrossprod(rep(1,n), h_mu)

  # Constant (but j-specific) AR(1) coef
  h_phi_all = tcrossprod(rep(1,n), h_phi)

  # Linear term:
  linht = matrix((ystar - m_st_all - h_mu_all)/v_st2_all)

  # Evolution precision matrix (n x p)
  evol_prec_mat = matrix(0, nrow = n, ncol = p);
  evol_prec_mat[1,] = 1/h_sigma_eta_0^2;
  evol_prec_mat[-1,] = 1/h_sigma_eta_t^2;

  # Lagged version, with zeros as appropriate (needed below)
  evol_prec_lag_mat = matrix(0, nrow = n, ncol = p);
  evol_prec_lag_mat[1:(n-1),] = evol_prec_mat[-1,]

  # Diagonal of quadratic term:
  Q_diag = matrix(1/v_st2_all +  evol_prec_mat + h_phi_all^2*evol_prec_lag_mat)

  # Off-diagonal of quadratic term:
  Q_off = matrix(-h_phi_all*evol_prec_lag_mat)[-(n*p)]

  # Quadratic term:
  QHt_Matrix = Matrix::bandSparse(n*p, k = c(0,1), diagonals = list(Q_diag, Q_off), symmetric = TRUE)

  # Cholesky:
  chQht_Matrix = Matrix::chol(QHt_Matrix)

  # Sample the log-vols:
  hsamp = h_mu_all + matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(length(linht))), nrow = n)


  # Return the (uncentered) log-vols
  hsamp
}

#' Sample components from a discrete mixture of normals
#'
#' Sample Z from 1,2,...,k, with P(Z=i) proportional to q[i]N(mu[i],sig2[i]).
#'
#' @param y Single numeric observation
#' @param mu vector of component means
#' @param sig vector of component standard deviations
#' @param q vector of component weights
#' @return Sample from {1,...,k}
#----------------------------------------------------------------------------
ncind <- function(y,mu,sig,q){
  sample(1:length(q),
         size = 1,
         prob = q*dnorm(y,mu,sig))
}

#----------------------------------------------------------------------------
#' Sampler for the stochastic volatility parameters
#'
#' Compute one draw of the normal stochastic volatility parameters.
#' The model assumes an AR(1) for the log-volatility.
#'
#' @param omega \code{n x p} matrix of errors
#' @param svParams list of parameters to be updated
#' @return List of relevant components in \code{svParams}: \code{sigma_wt}, the \code{n x p} matrix of standard deviations,
#' and additional parameters associated with SV model.
#'
#' @import stochvol
#' @export
sampleSVparams <- function(omega, svParams){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  for(j in 1:p){
    # First, check for numerical issues:
    svInput = omega[,j]; #if(all(svInput==0)) {svInput = 10^-8} else svInput = svInput + sd(svInput)/10^8

    # Sample the SV parameters:
    svsamp = stochvol::svsample_fast_cpp(svInput,
                                         startpara = list(
                                           mu = svParams$svParams[1,j],
                                           phi = svParams$svParams[2,j],
                                           sigma = svParams$svParams[3,j]),
                                         startlatent = svParams$ht[,j])# ,priorphi = c(10^4, 10^4));
    # Update the parameters:
    svParams$svParams[,j] = svsamp$para[1:3];
    svParams$ht[,j] = svsamp$latent
  }
  # Finally, up the evolution error SD:
  svParams$sigma_wt = exp(svParams$ht/2)

  # Check for numerically large values:
  svParams$sigma_wt[which(svParams$sigma_wt > 10^3, arr.ind = TRUE)] = 10^3

  return(svParams)
}

#----------------------------------------------------------------------------
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Printed timing message when reporting is due; otherwise NULL.
computeTimeRemaining <- function(nsi, timer0, nsims, nrep=1000){

  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==20) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
}

# ============================================================================
# Univariate sampler with AR/DAR support
# ============================================================================

#' MCMC sampler for univariate Bayesian trend filtering
#'
#' Fit Gaussian observations with AR, first-difference (D1), differenced AR
#' (DAR), or second-difference (D2) state evolution. Evolution errors have
#' DHS, HS, BL, SV or NIG priors. The Gibbs sampler alternates state-path,
#' shrinkage and observation-variance draws, and imputes missing responses.
#'
#' @param y Numeric length-n response vector; NA values are imputed. At least
#'   two observed values with positive variance and four time points (five for
#'   D2) are required by the retained evolution-prior initialisation/updates.
#' @param evol_error Evolution-error prior: DHS, HS, BL, SV, or NIG.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param useObsSV Whether observation errors have stochastic volatility.
#' @param nsave Positive integer number of draws to retain; at least two for DIC.
#' @param nburn Non-negative integer number of burn-in draws.
#' @param nskip Non-negative integer number skipped between saved draws.
#' @param mcmc_params Character vector or list of names: mu, yhat,
#'   evol_sigma_t2, obs_sigma_t2, dhs_phi, dhs_mean, or rho.
#' @param computeDIC Whether to compute conditional DIC using observed responses.
#' @param verbose Whether to print progress and elapsed time.
#' @return Named list of requested draws. Paths have dimensions nsave by n;
#'   scalar parameters have length nsave. Inapplicable requested parameters
#'   remain NULL. loglike is always returned; DIC and p_d each contain two
#'   versions when requested. A run_info attribute records settings and versions.
#' @note AR/DAR use a Uniform(0, 0.99) rho prior. These are smoothing draws
#'   conditional on all supplied data. See .state_operator() for boundary priors.
btf <- function(
    y,
    evol_error = "DHS",
    differencing_option = "D1",
    useObsSV = FALSE,
    nsave = 1000,
    nburn = 1000,
    nskip = 4,
    mcmc_params = list(
      "mu", "yhat", "evol_sigma_t2", "obs_sigma_t2",
      "dhs_phi", "dhs_mean", "rho"
    ),
    computeDIC = TRUE,
    verbose = TRUE) {

  # Validate controls before allocating arrays or entering the MCMC.
  call <- match.call()
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  evol_error <- match.arg(toupper(evol_error), c("DHS", "HS", "BL", "SV", "NIG"))
  differencing_option <- match.arg(differencing_option, c("AR", "D1", "DAR", "D2"))
  mcmc_params <- .check_mcmc(y, differencing_option,
    list(nsave = nsave, nburn = nburn, nskip = nskip, useObsSV = useObsSV,
         computeDIC = computeDIC, verbose = verbose),
    mcmc_params, regression = FALSE)

  n <- length(y)
  # Interpolation supplies starting values only; MCMC re-imputes missing y.
  observed <- !is.na(y)
  t01 <- seq(0, 1, length.out = n)

  is.missing <- which(is.na(y))
  any.missing <- length(is.missing) > 0
  y <- approxfun(t01, y, rule = 2)(t01)

  sigma_e <- sd(y, na.rm = TRUE)
  sigma_et <- rep(sigma_e, n)

  if (differencing_option %in% c("AR", "DAR")) {
    rho <- sample_AR1_param(
      beta = NULL,
      differencing_option = differencing_option,
      prior_type = "uniform"
    )
  } else {
    rho <- NULL
  }

  chol0 <- initChol.spam(
    n = n,
    differencing_option = differencing_option,
    rho = rho
  )

  mu <- sampleBTF(
    y,
    obs_sigma_t2 = sigma_et^2,
    evol_sigma_t2 = 0.01 * sigma_et^2,
    differencing_option = differencing_option,
    rho = rho,
    chol0 = chol0
  )

  omega <- compute_omega(mu, differencing_option, rho)
  mu0 <- compute_beta0(mu, differencing_option)

  evolParams <- initEvolParams(omega, evol_error = evol_error)
  evolParams0 <- initEvol0(mu0)

  if (useObsSV) {
    svParams <- initSV(y - mu)
    sigma_et <- as.numeric(svParams$sigma_wt)
  }

  # Allocate requested draws; DIC additionally needs mean and variance draws.
  mcmc_output <- vector("list", length(mcmc_params))
  names(mcmc_output) <- mcmc_params

  if (!is.na(match("mu", mcmc_params)) || computeDIC) {
    post_mu <- array(NA, c(nsave, n))
  }
  if (!is.na(match("yhat", mcmc_params))) {
    post_yhat <- array(NA, c(nsave, n))
  }
  if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) {
    post_obs_sigma_t2 <- array(NA, c(nsave, n))
  }
  if (!is.na(match("evol_sigma_t2", mcmc_params))) {
    post_evol_sigma_t2 <- array(NA, c(nsave, n))
  }
  if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == "DHS") {
    post_dhs_phi <- numeric(nsave)
  }
  if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == "DHS") {
    post_dhs_mean <- numeric(nsave)
  }
  if (differencing_option %in% c("AR", "DAR")) {
    post_rho <- numeric(nsave)
  }
  post_loglike <- numeric(nsave)

  nstot <- nburn + (nskip + 1) * nsave
  skipcount <- 0
  isave <- 0

  if (verbose) timer0 <- proc.time()[3]

  # One Gibbs sweep: missing y, rho, coefficient paths, shrinkage, then noise.
  for (nsi in seq_len(nstot)) {
    if (any.missing) {
      y[is.missing] <- mu[is.missing] +
        sigma_et[is.missing] * rnorm(length(is.missing))
    }

    # Initial-state priors precede the evolution-error variances in H beta.
    initial_rows <- if (differencing_option == "D2") 2L else 1L
    evol_sigma_t2 <- c(
      rep(evolParams0$sigma_w0^2, length.out = initial_rows),
      evolParams$sigma_wt^2
    )

    if (differencing_option %in% c("AR", "DAR")) {
      rho <- sample_AR1_param(
        beta = mu,
        differencing_option = differencing_option,
        evol_sigma_t2 = evol_sigma_t2,
        rho_current = rho,
        prior_type = "uniform"
      )
    }

    mu <- sampleBTF(
      y,
      obs_sigma_t2 = sigma_et^2,
      evol_sigma_t2 = evol_sigma_t2,
      differencing_option = differencing_option,
      rho = rho,
      chol0 = chol0
    )

    omega <- compute_omega(mu, differencing_option, rho)
    mu0 <- compute_beta0(mu, differencing_option)
    evolParams0 <- sampleEvol0(mu0, evolParams0, A = 1)

    # Update shrinkage scales and either time-varying or constant noise.
    if (useObsSV) {
      evolParams <- sampleEvolParams(
        omega,
        evolParams,
        1 / sqrt(n),
        evol_error
      )
      svParams <- sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et <- as.numeric(svParams$sigma_wt)
    } else {
      evolParams <- sampleEvolParams(
        omega,
        evolParams,
        sigma_e / sqrt(n),
        evol_error
      )

      if (evol_error == "DHS") {
        sigma_e <- uni.slice(sigma_e, g = function(x) {
          -(n + 2) * log(x) -
            0.5 * sum((y - mu)^2, na.rm = TRUE) / x^2 -
            log(1 + (sqrt(n) * exp(evolParams$dhs_mean0 / 2) / x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if (evol_error == "HS") {
        sigma_e <- 1 / sqrt(rgamma(
          n = 1,
          shape = n / 2,
          rate = sum((y - mu)^2, na.rm = TRUE) / 2
        ))
      }
      if (evol_error == "BL") {
        sigma_e <- 1 / sqrt(rgamma(
          n = 1,
          shape = n / 2 + length(evolParams$tau_j) / 2,
          rate = sum((y - mu)^2, na.rm = TRUE) / 2 +
            n * sum((omega / evolParams$tau_j)^2) / 2
        ))
      }
      if (evol_error %in% c("NIG", "SV")) {
        sigma_e <- 1 / sqrt(rgamma(
          n = 1,
          shape = n / 2,
          rate = sum((y - mu)^2, na.rm = TRUE) / 2
        ))
      }

      sigma_et <- rep(sigma_e, n)
    }

    # Retain every (nskip + 1)-th draw after burn-in.
    if (nsi > nburn) {
      skipcount <- skipcount + 1

      if (skipcount > nskip) {
        isave <- isave + 1

        if (!is.na(match("mu", mcmc_params)) || computeDIC) {
          post_mu[isave, ] <- mu
        }
        if (!is.na(match("yhat", mcmc_params))) {
          post_yhat[isave, ] <- mu + sigma_et * rnorm(n)
        }
        if (!is.na(match("obs_sigma_t2", mcmc_params)) || computeDIC) {
          post_obs_sigma_t2[isave, ] <- sigma_et^2
        }
        if (!is.na(match("evol_sigma_t2", mcmc_params))) {
          post_evol_sigma_t2[isave, ] <- c(
            rep(evolParams0$sigma_w0^2, length.out = initial_rows),
            evolParams$sigma_wt^2
          )
        }
        if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == "DHS") {
          post_dhs_phi[isave] <- evolParams$dhs_phi
        }
        if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == "DHS") {
          post_dhs_mean[isave] <- evolParams$dhs_mean
        }
        if (differencing_option %in% c("AR", "DAR")) {
          post_rho[isave] <- rho
        }

        # Score observed responses only; imputations are latent variables.
        post_loglike[isave] <- sum(dnorm(
          y[observed],
          mean = as.numeric(mu)[observed],
          sd = sigma_et[observed],
          log = TRUE
        ))
        skipcount <- 0
      }
    }

    if (verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if (!is.na(match("mu", mcmc_params))) mcmc_output$mu <- post_mu
  if (!is.na(match("yhat", mcmc_params))) mcmc_output$yhat <- post_yhat
  if (!is.na(match("obs_sigma_t2", mcmc_params))) {
    mcmc_output$obs_sigma_t2 <- post_obs_sigma_t2
  }
  if (!is.na(match("evol_sigma_t2", mcmc_params))) {
    mcmc_output$evol_sigma_t2 <- post_evol_sigma_t2
  }
  if (!is.na(match("dhs_phi", mcmc_params)) && evol_error == "DHS") {
    mcmc_output$dhs_phi <- post_dhs_phi
  }
  if (!is.na(match("dhs_mean", mcmc_params)) && evol_error == "DHS") {
    mcmc_output$dhs_mean <- post_dhs_mean
  }
  if (differencing_option %in% c("AR", "DAR") &&
      !is.na(match("rho", mcmc_params))) {
    mcmc_output$rho <- post_rho
  }

  mcmc_output$loglike <- post_loglike

  if (computeDIC) {
    # Conditional DIC follows dsp's two p_d definitions, using observed y only.
    loglike_hat <- sum(dnorm(
      y[observed],
      mean = colMeans(post_mu)[observed],
      sd = colMeans(sqrt(post_obs_sigma_t2))[observed],
      log = TRUE
    ))
    p_d <- c(
      2 * (loglike_hat - mean(post_loglike)),
      2 * var(post_loglike)
    )
    mcmc_output$DIC <- -2 * loglike_hat + 2 * p_d
    mcmc_output$p_d <- p_d
  }

  if (verbose) {
    print(paste("Total time: ", round(proc.time()[3] - timer0), "seconds"))
  }

  # Record reproducibility information alongside, rather than among, draws.
  attr(mcmc_output, "run_info") <- list(
    call = call, n = n, observed_n = sum(observed),
    differencing_option = differencing_option, evol_error = evol_error,
    nsave = nsave, nburn = nburn, nskip = nskip,
    useObsSV = useObsSV, computeDIC = computeDIC, mcmc_params = mcmc_params,
    rng_kind = RNGkind(), rng_state_at_entry = rng_state,
    R_version = R.version.string,
    package_versions = vapply(c("Matrix", "spam", "stochvol", "BayesLogit", "truncdist"),
                              function(pkg) as.character(utils::packageVersion(pkg)), character(1)))
  mcmc_output
}

#----------------------------------------------------------------------------
#' Initialise a sparse Cholesky structure for univariate trend filtering
#' @param n Number of time points.
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar autocorrelation for AR/DAR.
#' @return A spam.chol.NgPeyton factor for reuse in sampleBTF().
#' @note Numeric factorisation is updated each draw; a changed sparsity pattern
#' is handled by the sampling helper. This function consumes no random numbers.
initChol.spam <- function(n, differencing_option = "D1", rho = NULL) {
  Q <- build_Q(rep(1, n), rep(1, n), differencing_option, rho)
  spam::chol.spam(.as_spam(.scale_precision(Q)$Q))
}

#----------------------------------------------------------------------------
#' Draw one univariate coefficient path under AR, D1, DAR, or D2 evolution
#' @param y Finite numeric length-n response vector; no missing values.
#' @param obs_sigma_t2 Positive finite length-n observation variances.
#' @param evol_sigma_t2 Positive finite length-n evolution variances, including
#'   the initial-state variance (two initial-state variances for D2).
#' @param differencing_option One of AR, D1, DAR, or D2.
#' @param rho Scalar autocorrelation for AR/DAR.
#' @param chol0 Optional spam structure from initChol.spam(); NULL uses Matrix.
#' @return Numeric n by 1 matrix of sampled states, preserving the original API.
sampleBTF <- function(y, obs_sigma_t2, evol_sigma_t2,
                      differencing_option = "D1", rho = NULL, chol0 = NULL) {
  .check_data(y)
  .check_variances(obs_sigma_t2, "obs_sigma_t2")
  if (length(obs_sigma_t2) != length(y) || length(evol_sigma_t2) != length(y)) stop("Variance vectors must have length(y) entries.")
  Q <- build_Q(obs_sigma_t2, evol_sigma_t2, differencing_option, rho)
  matrix(.sample_gaussian_precision(Q, y / as.numeric(obs_sigma_t2), chol0), ncol = 1)
}
