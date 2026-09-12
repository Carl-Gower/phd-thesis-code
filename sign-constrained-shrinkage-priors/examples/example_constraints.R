# Small synthetic examples for the installed R overlays.

this_file <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
}
component_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
source(file.path(component_dir, "use_overlays.R"))

set.seed(20260909)
n <- 50L
X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
y <- as.numeric(X %*% c(0.8, -0.7) + rnorm(n, sd = 0.5))
dat <- data.frame(y = y, x1 = X[, 1], x2 = X[, 2])

fit_bayesreg <- bayesreg::bayesreg(
  y ~ x1 + x2,
  data = dat,
  prior = "ridge",
  n.samples = 100,
  burnin = 50,
  thin = 1,
  n.cores = 1,
  pos_con = 1,
  neg_con = 2
)

fit_dl <- R2D2::dl(
  X, y,
  mcmc.n = 100,
  thin = 1,
  print = FALSE,
  pos_con = 1,
  neg_con = 2
)

fit_spike_slab <- BoomSpikeSlab::lm.spike(
  y ~ x1 + x2,
  niter = 100,
  data = dat,
  ping = 0,
  seed = 20260909,
  model.options = BoomSpikeSlab::SsvsOptions(adaptive.cutoff = Inf),
  pos_con = 2,
  neg_con = 3
)

stopifnot(
  all(fit_bayesreg$beta[1, ] >= 0),
  all(fit_bayesreg$beta[2, ] <= 0),
  all(fit_dl$beta[, 1] >= 0),
  all(fit_dl$beta[, 2] <= 0),
  all(fit_spike_slab$beta[, 2] >= 0),
  all(fit_spike_slab$beta[, 3] <= 0)
)

cat("All three examples respected the requested signs.\n")
