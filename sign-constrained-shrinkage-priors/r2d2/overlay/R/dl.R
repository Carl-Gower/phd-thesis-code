# Derived from R2D2 1.0 (commit e734639929abb60e616c114ac7fe4e2beb5c7f9d).
# Sign-constraint and numerical-stability modifications by Carl Gower, 2025.
# SPDX-License-Identifier: MIT

#' Gibbs Sampler for Dirichlet-Laplace Prior
#'
#' \code{dl} aims to generate the posterior samples for Dirichlet-Laplace
#' prior.
#'
#' @param x An n-by-p input matrix, each row of which is an observation vector.
#' @param y An n-by-one vector, representing the response variable.
#' @param hyper The values of hyperparameters \code{(a1,b1)} in the prior of
#'   sigma^2.
#' @param a.prior The value of hyperparameter in the prior, which can be a value
#'   between (1/max(n,p), 1/2). By default, it is set at 1/max(n,p).
#' @param mcmc.n Number of MCMC samples, by default, 10000.
#' @param thin Thinning parameter of the chain. The default is 1, representing
#'   no thinning.
#' @param eps Tolerance of convergence, by default, 1e-7.
#' @param print Boolean variable determining whether to print the progress of
#'   the MCMC iterations, i.e., which iteration the function is currently on.
#'   The default is TRUE.
#' @param pos_con Vector of indices of coefficients constrained to be positive.
#' @param neg_con Vector of indices of coefficients constrained to be negative.
#'   Indices in both arguments refer to columns of \code{x}.
#'
#' @return A list containing posterior draws for code{beta}, code{psi},
#'   code{phi}, code{tau}, and code{sigma2}.
#'
#' @export
dl <- function(x, y, hyper, a.prior, mcmc.n = 10000, thin = 1, eps = 1e-07, print = TRUE,
               pos_con = NULL, neg_con = NULL) {

  EPS = 1e+20
  # n, p, max(n,p)
  p <- ncol(x)
  n <- nrow(x)
  max.np <- max(n, p)

  # Check the optional coefficient constraints.
  check_constraint_indices <- function(index, argument) {
    if (is.null(index) || length(index) == 0L) return(NULL)
    if (!is.numeric(index) || anyNA(index) || any(index != as.integer(index))) {
      stop(argument, " must contain whole-number coefficient indices.")
    }
    index <- sort(unique(as.integer(index)))
    if (any(index < 1L | index > p)) {
      stop(argument, " indices must be between 1 and ", p, ".")
    }
    index
  }
  pos_con <- check_constraint_indices(pos_con, "pos_con")
  neg_con <- check_constraint_indices(neg_con, "neg_con")
  if (length(intersect(pos_con, neg_con)) > 0L) {
    stop("pos_con and neg_con must not overlap.")
  }

  # priors IG(a1,b1) for sigma2.
  if (missing(hyper)) {
    a1 <- 0.001
    b1 <- 0.001
  } else {
    a1 <- hyper$a1
    b1 <- hyper$b1
  }

  # the discrete uniform value for a.
  if (missing(a.prior)) {
    a.prior <-  1/max.np
  }

  # Define variables to store posterior samples
  beta = phi = psi = matrix(0, nrow = mcmc.n, ncol = p)
  tau = sigma2 = rep(0, mcmc.n)

  #------------------------------------------
  # Initial values.
  tem <- stats::coef(stats::lm(y ~ x - 1))
  tem[which(is.na(tem))] <- eps
  if (!is.null(pos_con)) tem[pos_con] <- pmax(tem[pos_con], 0)
  if (!is.null(neg_con)) tem[neg_con] <- pmin(tem[neg_con], 0)
  beta[1, ] <- tem
  sigma2[1] <- MCMCpack::rinvgamma(1, shape = a1 + n/2, scale = b1 + 0.5 * crossprod(y - x %*% beta[1, ]))

  phi[1, ] <- rep(1/p, p)
  psi[1, ] <- stats::rexp(p, 0.5)
  tau[1] <- stats::rgamma(1, shape = p * a.prior, rate = 0.5)

  XTX <- crossprod(x)  #t(X)%*%X
  XTY <- crossprod(x, y)  #t(X)%*%y

  # Now MCMC runnings!
  for (k in 2:mcmc.n) {
    for (j in 1:thin) {

      Z <- stats::rnorm(p, mean = 0, sd = 1)

      #------------------------------------------
      # (ii) Sample beta| psi, phi, tau, y, sigma2
      d.te <- 1/(psi[k - 1, ] * (phi[k - 1, ]^2) * (tau[k - 1]^2))
      ad.te <- abs(d.te)
      sd.te <- sign(d.te)
      inx.e <- which(ad.te < eps)
      inx.E <- which(ad.te > EPS)
      d.te[inx.e] <- eps * sd.te[inx.e]
      d.te[which(is.infinite((d.te)))] <- EPS
      d.te[inx.E] <- EPS * sd.te[inx.E]

      if (length(d.te) == 1) {
        Dinv <- d.te
      } else {
        Dinv <- diag(d.te)
      }

      Vinv <- XTX + Dinv

      # Retain the original block update when no signs are constrained.
      if (is.null(pos_con) && is.null(neg_con)) {
        # Original efficient sampler.
        A <- Vinv

        # ********************************#********************************
        temQ <- chol(Vinv, pivot = T, tol = 0)
        pivot <- attr(temQ, "pivot")
        temc <- temQ[, order(pivot)]
        b <- XTY + t(temc) %*% Z * sqrt(sigma2[k - 1])

        beta[k, ] <- solve(A, b, tol = 0)

      } else {
        # Coordinate-wise Gibbs update under the sign restrictions.

        # Precision matrix and linear term. Guard against an exact zero scale.
        sig2_safe <- max(sigma2[k - 1], 1e-16)
        Omega <- Vinv / sig2_safe
        b_vec <- XTY / sig2_safe

        # Start the sweep from the preceding draw.
        beta_current <- beta[k - 1, ]

        # Update each coefficient from its conditional distribution.
        for (j in 1:p) {
          # Conditional variance
          cond_var <- 1 / Omega[j, j]

          # Conditional mean: (b_j - sum_{i != j} Omega_ji * beta_i) / Omega_jj
          sum_cross <- sum(Omega[j, -j] * beta_current[-j])
          cond_mean <- (b_vec[j] - sum_cross) * cond_var
          cond_sd <- sqrt(max(cond_var, 0)) # Max acts as a safety for floating-point 0

          # Apply a one-sided bound where requested.
          lb_j <- -Inf
          ub_j <- Inf
          if (!is.null(pos_con) && (j %in% pos_con)) lb_j <- 0
          if (!is.null(neg_con) && (j %in% neg_con)) ub_j <- 0

          # Draw from the appropriate univariate truncated normal.
          if (cond_sd < 1e-12) {
            # A nearly degenerate conditional is placed at its valid limit.
            beta_current[j] <- min(max(cond_mean, lb_j), ub_j)
          } else {
            beta_current[j] <- truncnorm::rtruncnorm(1, a = lb_j, b = ub_j,
                                                     mean = cond_mean, sd = cond_sd)
          }
        }

        beta[k, ] <- beta_current
      }


      beta[k, which(is.na(beta[k, ]))] <- 0

      #------------------------------------------
      # (i) Sample sigma2| beta, y
      sigma2[k] <- MCMCpack::rinvgamma(1, shape = a1 + n/2 + p/2, scale = b1 + 0.5 * crossprod(y - x %*% beta[k, ]) + sum(beta[k,
      ] * d.te * beta[k, ])/2)

      sigma.k <- sqrt(sigma2[k])

      #------------------------------------------
      # (iii) Sample psi|phi,tau,beta
      mu.te <- (phi[k - 1, ]/abs(beta[k, ]) * sigma.k) * tau[k - 1]
      amu.te <- abs(mu.te)
      smu.te <- sign(mu.te)
      inx.e <- which(amu.te < eps)
      inx.E <- which(amu.te > EPS)
      mu.te[inx.E] <- EPS * smu.te[inx.E]
      mu.te[inx.e] <- eps * smu.te[inx.e]

      psi[k, ] <- 1/rinvgauss_new(p, mean = mu.te, shape = 1)

      #------------------------------------------
      # (iv) Sample tau|phi,beta
      chi.te <- 2 * sum(abs(beta[k, ])/phi[k - 1, ])/sigma.k
      achi.te <- abs(chi.te)
      schi.te <- sign(chi.te)
      inx.e <- which(achi.te < eps)
      inx.E <- which(achi.te > EPS)
      chi.te[inx.e] <- eps * schi.te[inx.e]
      chi.te[inx.E] <- EPS * schi.te[inx.E]
      chi.te[which(chi.te == 0)] <- eps
      tau[k] <- rgig1(n = 1, lambda = p * a.prior - p, chi = chi.te, psi = 1)

      #------------------------------------------
      # (v) Sample phi|beta
      TT <- rep(0, p)
      tem <- 2 * abs(beta[k, ])/sigma.k
      atem <- abs(tem)
      stem <- sign(tem)
      inx.e <- which(atem < eps)
      inx.E <- which(atem > EPS)
      tem[inx.e] <- eps * stem[inx.e]
      tem[inx.E] <- EPS * stem[inx.E]
      tem[which(tem == 0)] <- eps
      TT <- apply(matrix(tem, ncol = 1), 1, function(xx) return(rgig1(n = 1, lambda = a.prior - 1, chi = xx, psi = 1)))

      Ts <- sum(TT)
      phi[k, ] <- TT/Ts
    }

    if (print) {
      if (k%%500 == 0) {
        print(paste(c("The ", k, "th sample."), collapse = ""))
      }
    }

    if (sqrt(sum((beta[k, ] - beta[k - 1, ])^2)) < eps & sqrt(sum((psi[k, ] - psi[k - 1, ])^2)) < eps & sqrt(sum((phi[k, ] -
                                                                                                                  phi[k - 1, ])^2)) < eps & abs(tau[k] - tau[k - 1]) < eps & abs(sigma2[k] - sigma2[k - 1]) < eps)
      break
  }

  return(list(beta = beta[(1):mcmc.n, ], psi = psi[(1):mcmc.n, ], phi = phi[(1):mcmc.n, ], tau = tau[(1):mcmc.n], sigma2 = sigma2[(1):mcmc.n]))
}


rgig1 <- function(n = 1, lambda, chi, psi) {
  n = as.integer(n)
  lambda = as.double(lambda)
  chi = as.double(chi)
  psi = as.double(psi)
  GIGrvg::rgig(n, lambda, chi, psi)
}
