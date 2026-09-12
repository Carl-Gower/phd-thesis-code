# Select this repository's installed packages for the current R session.
local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (!length(ofiles)) stop("Run this helper with source('use_overlays.R').")
  component_dir <- dirname(normalizePath(tail(ofiles, 1L)[[1L]], mustWork = TRUE))
  default_library <- file.path(component_dir, ".r-library", R.version$platform,
                              paste(R.version$major,
                                    strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L],
                                    sep = "."))
  library_dir <- Sys.getenv("PHD_R_LIBRARY", unset = default_library)
  library_dir <- normalizePath(library_dir, mustWork = TRUE)
  packages <- c("bayesreg", "R2D2", "Boom", "BoomSpikeSlab")
  for (package in packages) {
    description <- file.path(library_dir, package, "DESCRIPTION")
    if (!file.exists(description) ||
        !identical(unname(read.dcf(description,
          fields = "Config/phd-thesis-code/overlay")[1L, 1L]), "sign-constraints-1")) {
      stop("Run install_overlays.R first. Missing modified package: ", package,
           call. = FALSE)
    }
    if (package %in% loadedNamespaces()) {
      loaded_path <- getNamespaceInfo(asNamespace(package), "path")
      if (!identical(normalizePath(loaded_path),
                     normalizePath(file.path(library_dir, package)))) {
        stop(package, " is already loaded from a different library. ",
             "Restart R once, then source use_overlays.R before loading packages. ",
             "There is no need to uninstall anything.", call. = FALSE)
      }
    }
  }
  .libPaths(c(library_dir, .libPaths()))
  message("Using the sign-constraint packages in: ", library_dir)
})
