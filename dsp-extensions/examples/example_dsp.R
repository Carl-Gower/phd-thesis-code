# Small reproducible example of the DAR regression sampler.
# Run from the dsp-extensions folder with: Rscript examples/example_dsp.R

this_file <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  file_arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
  if (length(file_arg)) file_arg[1] else "examples/example_dsp.R"
}
component_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
source(file.path(component_dir, "R", "DSP_rho.R"))

set.seed(20260908)
n <- 60L
x <- rnorm(n)
X <- cbind(intercept = 1, predictor = x)

# Generate two smoothly changing coefficient paths and a noisy response.
time <- seq(0, 1, length.out = n)
beta_true <- cbind(
  intercept = 0.3 + 0.2 * sin(2 * pi * time),
  predictor = 0.8 - 0.5 * time
)
y <- rowSums(X * beta_true) + rnorm(n, sd = 0.3)

fit <- btf_reg(
  y,
  X,
  evol_error = "DHS",
  differencing_option = "DAR",
  rho_structure = "variable_specific",
  use_backfitting = TRUE,
  nsave = 50,
  nburn = 50,
  nskip = 0,
  mcmc_params = c("beta", "rho"),
  computeDIC = FALSE,
  verbose = FALSE
)

beta_median <- apply(fit$beta, c(2, 3), median)
rho_median <- apply(fit$rho, 2, median)

stopifnot(
  identical(dim(beta_median), c(n, ncol(X))),
  length(rho_median) == ncol(X),
  all(is.finite(beta_median)),
  all(is.finite(rho_median))
)

cat("Saved coefficient draws:", dim(fit$beta)[1], "x", dim(fit$beta)[2], "x", dim(fit$beta)[3], "\n")
cat("Posterior median rho:", paste(round(rho_median, 3), collapse = ", "), "\n")
cat("Synthetic example completed successfully.\n")
