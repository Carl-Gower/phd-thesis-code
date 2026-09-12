# Derived from R2D2 1.0 (commit e734639929abb60e616c114ac7fe4e2beb5c7f9d).
# Sign-constraint and numerical-stability modifications by Carl Gower, 2025.
# SPDX-License-Identifier: MIT

#' Gibbs Sampler for Marginal R2-D2
#'
#' \code{r2d2marg} adopts the Gibbs sampling algorithm, aiming to obtain a
#' sequence of posterior samples, which are approximately from the Marginal
#' R2-D2 prior.
#'
#' @param x An n-by-p input matrix, each row of which is an observation vector.
#' @param y An n-by-one vector, representing the response variable.
#' @param hyper The values of the hyperparameters in the prior, i.e., the values
#'   of \code{(a_pi, b, a1, b1)}.
#'   \itemize{
#'   \item{\code{a_pi} is the concentration parameter of the Dirichlet
#'   distribution, which controls the local shrinkage of phi. Smaller values of
#'   a_pi lead to most phi close to zero; while larger values of a_pi lead to a
#'   more uniform phi.}
#'   \item{\code{b} is the second shape parameter of beta prior for the
#'   R-squared.}
#'   \item{\code{a1} and \code{b1} are shape and scale parameters of Inverse
#'   Gamma distribution on sigma^2.} }
#'   \code{hyper} is set to be (1/(p^(b/2)n^(b/2)logn), 0.5, 0.001, 0.001) by
#'   default.
#' @param mcmc.n Number of MCMC samples, by default, 10000.
#' @param eps Tolerance of convergence, by default, 1e-7.
#' @param thin Thinning parameter of the chain. The default is 1, representing
#'   no thinning.
#' @param print Boolean variable determining whether to print the progress of
#'   the MCMC iterations, i.e., which iteration the function is currently on.
#'   The default is TRUE.
#' @param pos_con Vector of indices of coefficients constrained to be positive.
#' @param neg_con Vector of indices of coefficients constrained to be negative.
#'   Indices in both arguments refer to columns of \code{x}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{beta: Matrix (mcmc.n * p) of posterior samples for beta.}
#'   \item{sigma2: Vector (mcmc.n) of posterior samples for sigma^2.}
#'   \item{psi: Matrix (mcmc.n * p) of posterior samples for psi.}
#'   \item{w: Vector (mcmc.n) of posterior samples for the total prior
#'   probability w.}
#'   \item{xi: Vector (mcmc.n) of posterior samples for xi.}
#' }
#'
#' @export
r2d2marg <- function(x, y, hyper, mcmc.n = 10000, eps = 1e-07, thin = 1, print = TRUE,
                     pos_con = NULL, neg_con = NULL) {

  #------------------------------------------
  EPS = 1e+20
  # 1. n, p, max(n,p)
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

  # 2. hyperparameter values (a_pi, b, a1, b1)
  if (missing(hyper)) {
    b <- 0.5
    a_pi <- 1/(p^(b/2) * n^(b/2) * log(n))
    a1 <- 0.001
    b1 <- 0.001
  } else {
    a_pi <- hyper$a_pi
    b <- hyper$b
    a1 <- hyper$a1
    b1 <- hyper$b1
  }

  a <- a_pi * p

  # 3. Define variables to store posterior samples
  keep.beta = keep.phi = keep.psi = matrix(0, nrow = mcmc.n, ncol = p)
  keep.w = keep.xi = keep.sigma2 = rep(0, mcmc.n)

  # 4. Initial values.
  ## beta
  tem <- stats::coef(stats::lm(y ~ x - 1))
  tem[which(is.na(tem))] <- eps
  if (!is.null(pos_con)) tem[pos_con] <- pmax(tem[pos_con], 0)
  if (!is.null(neg_con)) tem[neg_con] <- pmin(tem[neg_con], 0)
  beta <- tem
  ## sigma
  sigma2 <- MCMCpack::rinvgamma(1, shape = a1 + n/2, scale = b1 + 0.5 * crossprod(y - x %*% beta))
  sigma <- sqrt(sigma2)
  ## phi
  phi <- rep(a_pi, p)
  ## xi
  xi <- a/(2 * b)
  ## w
  w <- b
  ## psi
  psi <- stats::rexp(p, rate = 0.5)
  #------------------------------------------

  XTX <- crossprod(x)  #t(X)%*%X
  XTY <- crossprod(x, y)  #t(X)%*%y

  keep.beta[1, ] <- beta
  keep.sigma2[1] <- sigma2
  keep.phi[1, ] <- phi
  keep.w[1] <- w
  keep.xi[1] <- xi
  keep.psi[1, ] <- psi

  #------------------------------------------
  #------------------------------------------
  # Now MCMC runnings!
  for (i in 2:mcmc.n) {
    for (j in 1:thin) {
      a <- a_pi * p

      # draw a number from N(0,1)
      Z <- stats::rnorm(p, mean = 0, sd = 1)

      #------------------------------------------
      # (i) Sample beta| phi, w, sigma2, y
      d.inv <- 1/(psi * phi * w)
      ad.inv <- abs(d.inv)
      sd.inv <- sign(d.inv)

      inx.e <- which(ad.inv < eps)
      inx.E <- which(ad.inv > EPS)
      d.inv[inx.e] <- eps * sd.inv[inx.e]
      d.inv[which(is.infinite((d.inv)))] <- EPS
      d.inv[inx.E] <- EPS * sd.inv[inx.E]

      if (length(d.inv) == 1) {
        Dinv <- d.inv
      } else {
        Dinv <- diag(d.inv)
      }

      Vinv <- XTX + Dinv

      # Retain the original block update when no signs are constrained.
      if (is.null(pos_con) && is.null(neg_con)) {
        # Original efficient sampler.

        # ********************************#********************************
        temQ <- chol(Vinv, pivot = T, tol = 1e-100000)
        pivot <- attr(temQ, "pivot")
        temc <- temQ[, order(pivot)]
        B <- XTY + t(temc) %*% Z * sigma

        # s = svd(Vinv) tes = s$u %*% diag(sqrt(s$d)) %*% t(s$v) b = XTY + t(tes)%*%Z*sqrt(sigma2[k-1])
        # b = XTY + t(chol(Vinv))%*%Z*sqrt(sigma2[k-1])
        # ********************************#********************************

        beta <- solve(Vinv, B, tol = 1e-10000)

      } else {
        # Multivariate truncated-normal update.

        # Remove small numerical asymmetries before factorisation.
        Vinv <- (Vinv + t(Vinv)) / 2

        # Compute the unscaled covariance.
        Q <- try(chol(Vinv), silent = TRUE)
        if (inherits(Q, "try-error")) {
          # Fall back to a direct solve if the Cholesky factorisation fails.
          Sigma_unscaled <- solve(Vinv, tol = 1e-12)
        } else {
          Rinv <- backsolve(Q, diag(p))
          Sigma_unscaled <- Rinv %*% t(Rinv)   # Vinv^{-1}
        }

        Sigma_unscaled <- (Sigma_unscaled + t(Sigma_unscaled)) / 2

        mu_beta <- as.vector(Sigma_unscaled %*% XTY)

        # Scale the covariance by the current residual variance.
        Sigma_beta <- Sigma_unscaled * sigma2

        Sigma_beta <- (Sigma_beta + t(Sigma_beta)) / 2

        # Add the smallest necessary diagonal adjustment if rounding has made
        # the covariance non-positive-definite.
        eigvals <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
        if (any(eigvals <= 0)) {
          jitter <- abs(min(eigvals)) + 1e-8
          diag(Sigma_beta) <- diag(Sigma_beta) + jitter
        }

        # Set one-sided bounds for the requested coefficients.
        lb <- rep(-Inf, p)
        ub <- rep(Inf, p)
        if (!is.null(pos_con)) lb[pos_con] <- 0
        if (!is.null(neg_con)) ub[neg_con] <- 0

        beta_draw <- TruncatedNormal::rtmvnorm(
          n    = 1,
          mu   = mu_beta,
          sigma = Sigma_beta,
          lb   = lb,
          ub   = ub
        )

        beta <- as.numeric(beta_draw)
      }

      beta[which(is.na(beta))] <- 0

      #------------------------------------------
      # (ii) Sample sigma2| beta, phi, w, y
      sigma2 <- MCMCpack::rinvgamma(1,
                                    shape = a1 + n/2 + p/2,
                                    scale = b1 + 0.5 * crossprod(y - x %*% beta) + sum(beta * d.inv * beta)/2
      )

      if (!is.finite(sigma2) || sigma2 <= 0) {
        warning("Non-finite sigma2 in iteration ", i, "; reverting to previous draw.")
        sigma2 <- keep.sigma2[i - 1]
      }
      sigma <- sqrt(sigma2)


      #------------------------------------------
      # (iii) Sample psi | beta, phi, w, sigma2
      #------------------------------------------
      mu.te <- sigma * sqrt(phi * w) / abs(beta)

      ## clamp mu.te away from 0 and Inf, and fix any non-finite values
      mu.te[!is.finite(mu.te)] <- eps
      mu.te <- pmin(pmax(mu.te, eps), EPS)

      ## rinvgauss_new: self-defined function
      psi <- 1 / rinvgauss_new(p, mean = mu.te, shape = 1)

      ## guard against bad psi values
      if (!all(is.finite(psi)) || any(psi <= 0)) {
        bad <- which(!is.finite(psi) | psi <= 0)
        if (length(bad) == length(psi)) {
          ## everything is bad – fall back to something sane
          psi[] <- 1
        } else {
          good_val <- min(psi[is.finite(psi) & psi > 0], na.rm = TRUE)
          psi[bad] <- good_val
        }
      }


      #------------------------------------------
      # (iv) Sample w | beta, phi, xi, sigma2
      #------------------------------------------

      ## make sure phi is positive and normalised
      if (!all(is.finite(phi)) || any(phi <= 0)) {
        bad_phi <- which(!is.finite(phi) | phi <= 0)
        if (length(bad_phi) == length(phi)) {
          ## all phi are bad – fall back to previous phi if available
          if (i > 2) {
            phi <- keep.phi[i - 1, ]
          } else {
            phi[] <- 1 / p
          }
        } else {
          phi[bad_phi] <- min(phi[is.finite(phi) & phi > 0], na.rm = TRUE)
        }
      }
      phi <- phi / sum(phi)

      ## compute chi.te safely
      chi.te <- sum(beta^2 / (psi * phi)) / sigma2

      if (!is.finite(chi.te) || chi.te <= 0) {
        chi.te <- eps
      }

      ## original bounding logic (now acting on a guaranteed-positive chi.te)
      achi.te <- abs(chi.te)
      schi.te <- sign(chi.te)
      inx.e <- which(achi.te < eps)
      inx.E <- which(achi.te > EPS)
      chi.te[inx.e] <- eps * schi.te[inx.e]
      chi.te[inx.E] <- EPS * schi.te[inx.E]
      chi.te[chi.te == 0] <- eps

      w <- GIGrvg::rgig(n = 1, lambda = a - p/2, chi = chi.te, psi = 4 * xi)



      #------------------------------------------
      # (v) Sample xi|w
      xi <- stats::rgamma(1, shape = a + b, rate = 1 + 2 * w)



      #------------------------------------------
      # (vi) Sample phi|beta, xi, sigma2, y
      TT <- rep(0, p)
      tem <- (beta^2 / sigma2) / psi

      ## make sure chi argument for GIG is valid
      tem[!is.finite(tem) | tem <= 0] <- eps

      atem <- abs(tem)
      stem <- sign(tem)
      inx.e <- which(atem < eps)
      inx.E <- which(atem > EPS)
      tem[inx.e] <- eps * stem[inx.e]
      tem[inx.E] <- EPS * stem[inx.E]
      tem[tem == 0] <- eps

      TT <- apply(matrix(tem, ncol = 1), 1, function(xx) {
        GIGrvg::rgig(n = 1, lambda = a_pi - 0.5, chi = xx, psi = 4 * xi)
      })

      Ts <- sum(TT)
      phi <- TT/Ts

    }

    if (print) {
      if (i%%500 == 0) {
        print(paste(c("The ", i, "th sample."), collapse = ""))
      }
    }

    keep.beta[i, ] <- beta
    keep.sigma2[i] <- sigma2
    keep.phi[i, ] <- phi
    keep.w[i] <- w
    keep.xi[i] <- xi
    keep.psi[i, ] <- psi

    if (sqrt(sum((keep.beta[i, ] - keep.beta[i - 1, ])^2)) < eps & abs(keep.xi[i] - keep.xi[i - 1]) <
        eps & sqrt(sum((keep.phi[i, ] - keep.phi[i - 1, ])^2)) < eps & abs(keep.w[i] - keep.w[i -
                                                                                              1]) < eps & abs(keep.sigma2[i] - keep.sigma2[i - 1]) < eps & sqrt(sum((keep.psi[i, ] - keep.psi[i -
                                                                                                                                                                                              1, ])^2)) < eps)
      break

  }

  return(list(beta = keep.beta, sigma2 = keep.sigma2, psi = keep.psi, w = keep.w, xi = keep.xi))
}
