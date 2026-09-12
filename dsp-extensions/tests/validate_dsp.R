# Reproducible validation for R/DSP_rho.R
# Run from the component root: Rscript tests/validate_dsp.R
# An alternative source path may be supplied as the first argument.
# No additional test packages are required.
# Tests use synthetic data; they validate implementation consistency rather
# than MCMC convergence for an empirical application. A failed assertion stops.

args <- commandArgs(trailingOnly = TRUE)
this_file <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else {
  file_arg <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
  if (length(file_arg)) file_arg[1] else "DSP_validation.R"
}
script <- if (length(args)) args[1] else file.path(dirname(this_file), "..", "R", "DSP_rho.R")
source(script)
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAILED: ", label)
  checks <<- checks + 1L
}
close_to <- function(actual, expected, label, tolerance = 1e-9) {
  check(max(abs(as.numeric(actual) - as.numeric(expected))) <= tolerance, label)
}
expect_error <- function(expr, pattern) {
  result <- tryCatch({force(expr); NULL}, error = function(e) conditionMessage(e))
  check(!is.null(result) && grepl(pattern, result, fixed = TRUE), paste("reject", pattern))
}

# Construct an independent dense operator by applying the written residual
# equations to each basis vector. This deliberately does not call build_Q,
# compute_omega, .state_operator or any precision builder under test.
dense_H <- function(n, p, option, rho) {
  apply_equations <- function(v) {
    B <- matrix(v, n, p, byrow = TRUE)
    E <- B
    for (j in seq_len(p)) for (t in 2:n) {
      E[t,j] <- switch(option,
        AR = B[t,j] - rho[j]*B[t-1,j],
        D1 = B[t,j] - B[t-1,j],
        DAR = if (t == 2) B[t,j]-(1-rho[j])*B[t-1,j] else
          B[t,j]-(1+rho[j])*B[t-1,j]+rho[j]*B[t-2,j],
        D2 = if (t == 2) B[t,j] else B[t,j]-2*B[t-1,j]+B[t-2,j])
    }
    as.vector(t(E))
  }
  vapply(seq_len(n*p), function(k) apply_equations(diag(n*p)[,k]), numeric(n*p))
}
dense_posterior <- function(X, y, obs, evol, option, rho) {
  n <- nrow(X); p <- ncol(X)
  Z <- matrix(0, n, n*p)
  for(t in seq_len(n)) Z[t, ((t-1)*p+1):(t*p)] <- X[t,]
  H <- dense_H(n,p,option,rho)
  Q <- crossprod(Z, Z / obs) + crossprod(H, H / as.vector(t(evol)))
  b <- as.vector(crossprod(Z, y/obs))
  list(Q=Q, mean=as.vector(solve(Q,b)), covariance=solve(Q))
}

set.seed(20260907)
for (n in c(3L, 8L)) for (p in c(1L, 2L, 4L)) {
  X <- matrix(rnorm(n*p), n, p)
  X[1,] <- 0
  if (p > 1) X[,2] <- X[,1] # deliberately collinear; priors are proper
  y <- rnorm(n); obs <- runif(n,.4,2); evol <- matrix(runif(n*p,.05,.8),n)
  for(option in c("AR","D1","DAR","D2")) for (structure in c("shared","variable_specific")) {
    rho <- if (structure == "shared") rep(.55,p) else seq(-.7,.8,length.out=p)
    rho_arg <- if(option %in% c("AR","DAR")) if(structure == "shared") rho[1] else rho else NULL
    ref <- dense_posterior(X,y,obs,evol,option,rho)
    Q <- .reg_precision(obs,evol,build_XtX(X),option,rho_arg)
    close_to(Q, ref$Q, paste("dense precision",n,p,option,structure))
    check(min(eigen(ref$Q,symmetric=TRUE,only.values=TRUE)$values)>0, "proper precision is positive definite")
    H <- dense_H(n,p,option,rho)
    B <- matrix(rnorm(n*p),n,p)
    E <- matrix(H %*% as.vector(t(B)),n,p,byrow=TRUE)
    first <- if(option == "D2") 2 else 1
    close_to(compute_omega(B,option,rho_arg), E[-seq_len(first),,drop=FALSE], "residual equations")
    close_to(compute_beta0(B,option), B[seq_len(first),,drop=FALSE], "initial states")
    # With z=0, the implemented solve must return Q^{-1}b. With each unit
    # vector, its innovation transform must have covariance Q^{-1}.
    sampler <- .sample_gaussian_precision
    environment(sampler) <- new.env(parent=environment(.sample_gaussian_precision))
    environment(sampler)$rnorm <- function(n) rep(0,n)
    b <- as.vector(t(X * (y/obs)))
    centre <- sampler(Q,b)
    close_to(centre,ref$mean,"Gaussian conditional mean")
    innovation <- vapply(seq_len(n*p),function(k) {
      environment(sampler)$rnorm <- function(n) as.numeric(seq_len(n)==k)
      sampler(Q,b)-centre
    },numeric(n*p))
    close_to(tcrossprod(innovation),ref$covariance,"Gaussian conditional covariance",1e-8)
    for (j in seq_len(p)) {
      Qj <- build_Q(obs,evol[,j],option,if(is.null(rho_arg))NULL else rho[j])
      Hj <- dense_H(n,1,option,rho[j])
      close_to(Qj,diag(1/obs)+crossprod(Hj,Hj/evol[,j]),"univariate precision")
    }
  }
}
cat("PASS: dense model equations, boundary rows, means and covariances\n")

# Finite draws at the minimum low-level sizes, including DAR with only two states.
for(option in c("AR","D1","DAR","D2")) {
  n <- if(option == "D2") 3L else 2L
  rho <- if(option %in% c("AR","DAR")) .6 else NULL
  check(all(is.finite(sampleBTF(seq_len(n),rep(1,n),rep(.2,n),option,rho))),"short univariate path")
  X <- matrix(seq_len(n),n,1)
  check(all(is.finite(sampleBTF_reg(seq_len(n),X,rep(1,n),matrix(.2,n,1),build_XtX(X),option,rho))),"short joint path")
}

# Exercise the actual spam draw and its permutations statistically against an
# independent dense Gaussian. Tolerances are in Monte Carlo standard errors.
set.seed(719)
n <- 4L; p <- 2L; draws <- 2000L
X <- cbind(1,c(-.8,.3,1.1,-.2)); y <- c(.2,-.1,.6,.7)
obs <- c(.3,.8,.5,.4); evol <- matrix(c(.8,.4,.5,.7,.5,.3,.8,.6),n,p)
for(option in c("AR","D1","DAR","D2")) {
  rho <- if(option %in% c("AR","DAR")) c(.2,.8) else NULL
  ref <- dense_posterior(X,y,obs,evol,option,if(is.null(rho))rep(0,p)else rho)
  cache <- initCholReg.spam(obs,evol,build_XtX(X),option,rho)
  sims <- replicate(draws,as.vector(t(sampleBTF_reg(y,X,obs,evol,build_XtX(X),option,rho,cache))))
  mean_error <- abs(rowMeans(sims)-ref$mean)/sqrt(diag(ref$covariance)/draws)
  cov_se <- sqrt((ref$covariance^2+outer(diag(ref$covariance),diag(ref$covariance)))/(draws-1))
  check(max(mean_error)<6,"spam Monte Carlo means")
  check(max(abs(cov(t(sims))-ref$covariance)/cov_se)<6,"spam Monte Carlo covariance")
}
cat("PASS: spam joint draws against dense Gaussian moments (8,000 draws)\n")

# A random-order Gibbs sweep must preserve the joint conditional distribution.
# Start every replicate at an independent exact joint draw, then sweep once.
set.seed(151)
for(option in c("AR","D1","DAR","D2")) {
  rho <- if(option %in% c("AR","DAR")) c(.2,.8) else NULL
  ref <- dense_posterior(X,y,obs,evol,option,if(is.null(rho))rep(0,p)else rho)
  root_cov <- t(chol(ref$covariance))
  sims <- replicate(draws, {
    beta0 <- matrix(ref$mean+root_cov%*%rnorm(n*p),n,p,byrow=TRUE)
    as.vector(t(sampleBTF_reg_backfit(y,X,beta0,obs,evol,option,rho)))
  })
  check(max(abs(rowMeans(sims)-ref$mean)/sqrt(diag(ref$covariance)/draws))<6,"backfitting preserves means")
  cov_se <- sqrt((ref$covariance^2+outer(diag(ref$covariance),diag(ref$covariance)))/(draws-1))
  check(max(abs(cov(t(sims))-ref$covariance)/cov_se)<6,"backfitting preserves covariance")
}
cat("PASS: backfitting leaves the joint conditional invariant (8,000 sweeps)\n")

# Check rho's conditional CDF against numerical integration of the original
# residual likelihood, not against a duplicate weighted-regression formula.
set.seed(619)
B <- cbind(c(1,1.4,.8,1.8),c(-.8,-1.4,-.4,-1.3))
Es <- matrix(c(.5,.1,.2,.4,.7,.2,.3,.6),4,2)
for(option in c("AR","DAR")) for(prior in c("uniform","truncated_normal")) {
  log_density <- function(r) {
    E <- matrix(dense_H(4,2,option,rep(r,2))%*%as.vector(t(B)),4,2,byrow=TRUE)
    -sum(E[-1,]^2/Es[-1,])/2 + if(prior=="uniform") 0 else dnorm(r,.3,.7,log=TRUE)
  }
  grid <- seq(0,.99,length.out=201)
  offset <- max(vapply(grid,log_density,numeric(1)))
  density <- function(r) vapply(r,function(z) exp(log_density(z)-offset),numeric(1))
  target <- integrate(density,0,.5)$value/integrate(density,0,.99)$value
  r <- replicate(3000,sample_AR1_param(B,option,Es,prior_type=prior,mu_prior=.3,sigma_prior=.7))
  check(all(is.finite(r)&r>=0&r<=.99),"rho bounds")
  check(abs(mean(r<=.5)-target)<6*sqrt(target*(1-target)/length(r))+.003,"rho CDF against quadrature")
}
# This likelihood has untruncated mean 3: clamping it to .99 changes the CDF.
Btail <- c(1,3)
rtail <- replicate(5000,sample_AR1_param(Btail,"AR",c(1,.2)))
den <- function(r) exp(-((3-r)^2-(3-.99)^2)/.4)
target <- integrate(den,0,.9)$value/integrate(den,0,.99)$value
check(abs(mean(rtail<=.9)-target)<.035,"unclamped outside-boundary posterior")
check(all(is.finite(replicate(100,.rtrunc_normal(50,.1,0,.99)))),"extreme truncated-normal tail")
cat("PASS: AR/DAR rho conditionals, prior choices and tail behaviour\n")

# All complete-MCMC combinations, using short runs solely to test execution.
# Missing y, a zero predictor, an intercept and collinearity are included.
set.seed(932)
n <- 16L; X <- cbind(1,rnorm(n),0); y <- rnorm(n); y[c(1,8,16)] <- NA
X[,3] <- X[,2]; X[c(3,6),2:3] <- 0
fit_count <- 0L
for(option in c("AR","D1","DAR","D2")) for(prior in c("DHS","HS","BL","SV","NIG")) for(obs_sv in c(FALSE,TRUE)) {
  for(backfit in c(FALSE,TRUE)) for(structure in c("shared","variable_specific")) {
    fit <- btf_reg(y,X,evol_error=prior,differencing_option=option,rho_structure=structure,
                   useObsSV=obs_sv,nsave=3,nburn=2,nskip=1,use_backfitting=backfit,verbose=FALSE)
    check(identical(dim(fit$beta),c(3L,n,3L)),"regression saved dimensions")
    check(all(is.finite(unlist(fit))),"regression finite output")
    check(attr(fit,"run_info")$observed_n==13,"observed count metadata")
    mu <- fit$mu; sdobs <- sqrt(fit$obs_sigma_t2); observed <- !is.na(y)
    expected <- vapply(1:3,function(i)sum(dnorm(y[observed],mu[i,observed],sdobs[i,observed],log=TRUE)),numeric(1))
    close_to(fit$loglike,expected,"likelihood excludes imputed responses",1e-8)
    fit_count <- fit_count+1L
  }
  fit <- btf(y,evol_error=prior,differencing_option=option,useObsSV=obs_sv,
             nsave=3,nburn=2,nskip=1,verbose=FALSE)
  check(identical(dim(fit$mu),c(3L,n)),"univariate saved dimensions")
  check(all(is.finite(unlist(fit))),"univariate finite output")
  fit_count <- fit_count+1L
}
cat("PASS:",fit_count,"complete MCMC configurations across priors and observation SV\n")

# Regression initialisation must never call a joint sampler when backfitting.
env <- new.env(parent=globalenv())
for(nm in c("initCholReg.spam","sampleBTF_reg","build_XtX")) env[[nm]] <- function(...) stop("Joint path called")
backfit_wrapper <- btf_reg; environment(backfit_wrapper) <- env
set.seed(12)
check(all(is.finite(unlist(backfit_wrapper(rnorm(12),matrix(1,12,1),nsave=2,nburn=1,nskip=0,use_backfitting=TRUE,verbose=FALSE)))),"backfitting avoids joint initialisation")

# Output selection, wrapper routing, minimum supported MCMC sizes, RNG replay.
for(option in c("AR","D1","DAR","D2")) {
  n <- if(option=="D2")5L else 4L
  set.seed(154)
  yy <- rnorm(n)
  for(prior in c("DHS","HS","BL","SV","NIG")) {
    fit <- btf_reg(yy,matrix(1,n,1),evol_error=prior,differencing_option=option,
                   nsave=2,nburn=1,nskip=0,use_backfitting=TRUE,verbose=FALSE)
    check(all(is.finite(unlist(fit))),"minimum length one-predictor MCMC")
  }
}
set.seed(214)
yy <- rnorm(12); XX <- cbind(1,rnorm(12))
run_fit <- function() btf_reg(yy,XX,differencing_option="DAR",nsave=2,nburn=2,nskip=0,verbose=FALSE)
set.seed(982); first <- run_fit()
set.seed(982); second <- run_fit()
check(identical(first,second),"seeded reproducibility")
selected <- btf_reg(yy,X=NULL,nsave=1,nburn=0,nskip=0,mcmc_params=c("rho"),computeDIC=FALSE,verbose=FALSE)
check(identical(names(selected),c("rho","loglike")),"X=NULL output selection without beta")
expect_error(btf_reg(yy,XX[-1,],verbose=FALSE),"X must be")
expect_error(btf(yy,nsave=0),"nsave must be")
expect_error(btf(yy,nsave=1),"computeDIC requires")
expect_error(btf(rep(NA_real_,12)),"At least two observed")
expect_error(btf(rep(1,12)),"positive sample variance")
expect_error(btf(yy,evol_error="unknown"),"arg")
expect_error(btf(yy,mcmc_params="unknown"),"mcmc_params")
expect_error(build_Q(rep(1,4),c(1,0,1,1)),"evol_sigma_t2")
expect_error(build_Q(rep(1,4),rep(1,4),"AR",NULL),"rho must be")
expect_error(sample_AR1_param(yy,"AR",rep(1,11)),"matching n by p")
expect_error(.rho_vector(c(.1,.2),3,"AR"),"rho must be")
cat("PASS: validation, output selection, initialisation and reproducibility\n")
cat("ALL",checks,"CHECKS PASSED\n")
print(sessionInfo())
