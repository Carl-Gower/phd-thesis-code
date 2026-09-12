#!/usr/bin/env Rscript
# Download the pinned sources, apply the changes, and install in a local library.
local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  sourced <- length(ofiles) > 0L
  script_path <- if (sourced) tail(ofiles, 1L)[[1L]] else {
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
  }
  script_path <- normalizePath(script_path, mustWork = TRUE)
  component_dir <- dirname(script_path)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") {
    "Rscript.exe"
  } else "Rscript")

  # Installation always runs in a fresh process, including from RStudio.
  if (sourced) {
    status <- system2(rscript, c("--vanilla", shQuote(script_path)))
    if (status != 0L) stop("Installation failed; see the message above.", call. = FALSE)
  } else {
    if (getRversion() < "4.5.0") {
      stop("These versions of Boom and BoomSpikeSlab require R 4.5.0 or newer.",
           call. = FALSE)
    }
    r_series <- paste(R.version$major,
                      strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L], sep = ".")
    default_library <- file.path(component_dir, ".r-library", R.version$platform, r_series)
    library_dir <- Sys.getenv("PHD_R_LIBRARY", unset = default_library)
    dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
    library_dir <- normalizePath(library_dir, mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    # BoomSpikeSlab's Makevars invokes Rscript to locate Boom's headers and library.
    # Set R_LIBS as well as .libPaths so that subprocess finds the modified Boom.
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    Sys.setenv(PHD_R_LIBRARY = library_dir)
    options(timeout = max(600, getOption("timeout")))
    build_root <- Sys.getenv("PHD_R_BUILD", unset = file.path(component_dir, ".r-build"))
    dir.create(build_root, recursive = TRUE, showWarnings = FALSE)
    build_root <- normalizePath(build_root, mustWork = TRUE)
    download_dir <- file.path(build_root, "downloads")
    dir.create(download_dir, showWarnings = FALSE)
    build_dir <- tempfile("install-", tmpdir = build_root)
    dir.create(build_dir)
    cat("Installing into:", library_dir, "\nBuild logs:", build_dir, "\n")

    cran_urls <- function(package, version) {
      filename <- paste0(package, "_", version, ".tar.gz")
      c(paste0("https://cran.r-project.org/src/contrib/", filename),
        paste0("https://cran.r-project.org/src/contrib/Archive/", package, "/", filename))
    }
    sources <- list(
      bayesreg = list(file = "bayesreg_1.3.tar.gz", urls = cran_urls("bayesreg", "1.3"),
                      directory = "bayesreg", md5 = "bec9a66d3dc0d3b9da0dbe3680416510"),
      r2d2 = list(file = "R2D2_e7346399.tar.gz", urls = paste0(
        "https://codeload.github.com/yandorazhang/R2D2/tar.gz/",
        "e734639929abb60e616c114ac7fe4e2beb5c7f9d"),
        directory = "R2D2-e734639929abb60e616c114ac7fe4e2beb5c7f9d",
        md5 = "46c509163291b206ea317ab054960ffe"),
      boom = list(file = "Boom_0.9.16.tar.gz", urls = cran_urls("Boom", "0.9.16"),
                  directory = "Boom", md5 = "d3a299c3c467a7075d762562028b3907"),
      boomspikeslab = list(file = "BoomSpikeSlab_1.2.7.tar.gz",
        urls = cran_urls("BoomSpikeSlab", "1.2.7"), directory = "BoomSpikeSlab",
        md5 = "e97c76a79079403f8d01ccbfe465fa96")
    )
    source(file.path(component_dir, "apply_overlay.R"), local = TRUE)
    source_dirs <- character()
    for (component in names(sources)) {
      spec <- sources[[component]]
      archive <- file.path(download_dir, spec$file)
      valid_archive <- function() {
        file.exists(archive) && identical(unname(tools::md5sum(archive)), spec$md5)
      }
      if (!valid_archive()) {
        downloaded <- FALSE
        for (url in spec$urls) {
          status <- tryCatch(
            suppressWarnings(download.file(url, archive, mode = "wb", quiet = TRUE)),
            error = function(e) 1L)
          if (isTRUE(status == 0L) && valid_archive()) {
            downloaded <- TRUE
            break
          }
        }
        if (!downloaded) {
          stop("Could not download the verified source for ", component,
               ". Check your internet connection and the source URLs in this script.",
               call. = FALSE)
        }
      }
      untar(archive, exdir = build_dir)
      source_dir <- file.path(build_dir, spec$directory)
      apply_overlay(component, source_dir, component_dir)
      source_dirs[[component]] <- source_dir
    }

    dependencies <- c("pgdraw", "doParallel", "foreach", "TruncatedNormal", "MASS",
                      "statmod", "GIGrvg", "MCMCpack", "mvtnorm", "onlinePCA",
                      "chemometrics", "truncnorm")
    missing <- setdiff(dependencies, rownames(installed.packages()))
    if (length(missing)) {
      install.packages(missing, lib = library_dir, repos = "https://cloud.r-project.org")
      missing <- setdiff(dependencies, rownames(installed.packages()))
      if (length(missing)) stop("Dependencies did not install: ",
                                paste(missing, collapse = ", "), call. = FALSE)
    }
    r <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
    # This order is required: Boom must be installed before BoomSpikeSlab.
    for (component in names(source_dirs)) {
      log_file <- file.path(build_dir, paste0(component, "-install.log"))
      cat("Building", component, "...\n")
      status <- system2(r, c("CMD", "INSTALL", "--preclean",
                             paste0("--library=", shQuote(library_dir)),
                             shQuote(source_dirs[[component]])),
                        stdout = log_file, stderr = log_file)
      if (status != 0L) {
        cat(tail(readLines(log_file, warn = FALSE), 25L), sep = "\n")
        stop("Installation failed for ", component, ". Full log: ", log_file,
             call. = FALSE)
      }
    }
    test_file <- file.path(component_dir, "tests", "smoke_installed.R")
    status <- system2(rscript, c("--vanilla", shQuote(test_file)))
    if (status != 0L) stop("Installed-package checks failed.", call. = FALSE)
    cat("\nInstallation and checks passed. In R, run:\n")
    cat("source(", dQuote(file.path(component_dir, "use_overlays.R"), q = FALSE), ")\n",
        sep = "")
  }
})
