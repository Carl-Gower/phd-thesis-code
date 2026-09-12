# Smoke tests for the static sign-constraint overlays.

this_file <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  file_arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
  if (length(file_arg)) file_arg[1] else "tests/smoke_installed.R"
}
component_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)

source(file.path(component_dir, "use_overlays.R"))

required <- c("bayesreg", "R2D2", "Boom", "BoomSpikeSlab",
              "TruncatedNormal", "truncnorm")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing installed packages: ", paste(missing, collapse = ", "))

# Test the functions users actually call from the installed packages.
stopifnot(all(c("pos_con", "neg_con") %in% names(formals(bayesreg::bayesreg))),
          all(c("pos_con", "neg_con") %in% names(formals(R2D2::dl))),
          all(c("pos_con", "neg_con") %in% names(formals(R2D2::r2d2marg))),
          all(c("pos_con", "neg_con") %in% names(formals(BoomSpikeSlab::lm.spike))))
for (package in required[1:4]) {
  cat(package, as.character(packageVersion(package)), find.package(package), "\n")
}

expect_error <- function(expr, pattern) {
  message <- tryCatch({
    force(expr)
    NULL
  }, error = function(condition) conditionMessage(condition))
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

set.seed(20260909)
n <- 40L
X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
# Deliberately oppose the constraints to catch ignored sign arguments.
y <- as.numeric(X %*% c(-0.9, 0.8) + rnorm(n, sd = 0.5))
dat <- data.frame(y = y, x1 = X[, 1], x2 = X[, 2])

fit_bayesreg <- bayesreg::bayesreg(
  y ~ x1 + x2, dat,
  prior = "ridge", n.samples = 20, burnin = 10, thin = 1,
  n.cores = 1, pos_con = 1, neg_con = 2
)
stopifnot(all(fit_bayesreg$beta[1, ] >= 0),
          all(fit_bayesreg$beta[2, ] <= 0))

fit_dl <- R2D2::dl(X, y, mcmc.n = 20, print = FALSE,
                      pos_con = 1, neg_con = 2)
fit_r2d2 <- R2D2::r2d2marg(X, y, mcmc.n = 20, print = FALSE,
                              pos_con = 1, neg_con = 2)
stopifnot(all(fit_dl$beta[, 1] >= 0), all(fit_dl$beta[, 2] <= 0),
          all(fit_r2d2$beta[, 1] >= 0), all(fit_r2d2$beta[, 2] <= 0))

fit_spike_slab <- BoomSpikeSlab::lm.spike(
  y ~ x1 + x2, niter = 30, data = dat, ping = 0, seed = 20260909,
  model.options = BoomSpikeSlab::SsvsOptions(adaptive.cutoff = Inf),
  prior.inclusion.probabilities = rep(1, 3),
  pos_con = 2, neg_con = 3
)
stopifnot(all(fit_spike_slab$beta[, 2] > 0),
          all(fit_spike_slab$beta[, 3] < 0))

expect_error(
  R2D2::dl(X, y, mcmc.n = 3, print = FALSE,
              pos_con = 1, neg_con = 1),
  "must not overlap"
)
expect_error(
  bayesreg::bayesreg(y ~ x1 + x2, dat, n.samples = 2, burnin = 1,
                        thin = 1, n.cores = 1, pos_con = 3),
  "indices must be between"
)
expect_error(
  BoomSpikeSlab::lm.spike(
    y ~ x1 + x2, niter = 2, data = dat, ping = 0,
    model.options = BoomSpikeSlab::SsvsOptions(adaptive.cutoff = 1),
    pos_con = 2
  ),
  "classic SSVS"
)

cat("STATIC_CONSTRAINT_SMOKE_TESTS_PASSED\n")
