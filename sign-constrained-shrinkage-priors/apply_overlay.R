#!/usr/bin/env Rscript
# Apply modified files to an unpacked upstream source package.

apply_overlay <- function(component_name, source_directory, component_root) {
  components <- list(
    bayesreg = list(overlay = "bayesreg/overlay", package = "bayesreg",
                   version = "^1\\.3(\\.0)?$", patched = "1.3.0.9000"),
    r2d2 = list(overlay = "r2d2/overlay", package = "R2D2",
                version = "^1\\.0(\\.0)?$", patched = "1.0.0.9000"),
    boom = list(overlay = "boom-spike-slab/boom-0.9.16-overlay", package = "Boom",
                version = "^0\\.9\\.16$", patched = "0.9.16.9000"),
    boomspikeslab = list(
      overlay = "boom-spike-slab/boom-spike-slab-1.2.7-overlay",
      package = "BoomSpikeSlab", version = "^1\\.2\\.7$", patched = "1.2.7.9000")
  )
  component_name <- tolower(component_name)
  component <- components[[component_name]]
  if (is.null(component)) {
    stop("Unknown component: ", component_name, call. = FALSE)
  }
  source_directory <- normalizePath(source_directory, mustWork = TRUE)
  description_file <- file.path(source_directory, "DESCRIPTION")
  if (!file.exists(description_file) ||
      file.exists(file.path(source_directory, "Meta", "package.rds"))) {
    stop("Use an unpacked source package, not an installed package: ",
         source_directory, call. = FALSE)
  }
  description <- read.dcf(description_file)[1L, ]
  version <- unname(description[["Version"]])
  if (!identical(unname(description[["Package"]]), component$package) ||
      !(grepl(component$version, version) || identical(version, component$patched))) {
    stop("Expected ", component$package, " at the upstream version in README.md; found ",
         description[["Package"]], " ", version, ".", call. = FALSE)
  }
  overlay_directory <- normalizePath(file.path(component_root, component$overlay),
                                     mustWork = TRUE)
  overlay_files <- list.files(overlay_directory, recursive = TRUE, full.names = TRUE)
  if (!length(overlay_files)) stop("The overlay is empty.", call. = FALSE)
  relative_paths <- substring(overlay_files, nchar(overlay_directory) + 2L)
  target_files <- file.path(source_directory, relative_paths)
  if (!all(file.exists(target_files))) {
    stop("The source tree is missing: ",
         paste(relative_paths[!file.exists(target_files)], collapse = ", "),
         call. = FALSE)
  }
  copied <- file.copy(overlay_files, target_files, overwrite = TRUE)
  if (!all(copied)) {
    stop("Could not copy: ", paste(relative_paths[!copied], collapse = ", "),
         call. = FALSE)
  }
  # Distinguish the installed modifications from the upstream releases.
  metadata <- readLines(description_file, warn = FALSE)
  metadata <- sub("^Version:.*$", paste("Version:", component$patched), metadata)
  metadata <- metadata[!grepl("^Config/phd-thesis-code/overlay:", metadata)]
  writeLines(c(metadata, "Config/phd-thesis-code/overlay: sign-constraints-1"),
             description_file)
  cat("Applied", component_name, "overlay to", source_directory, "\n")
  cat(paste0("  ", relative_paths, "\n"), sep = "")
  invisible(source_directory)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop(paste("Usage: Rscript apply_overlay.R COMPONENT SOURCE_DIRECTORY",
               "Components: bayesreg, r2d2, boom, boomspikeslab", sep = "\n"),
         call. = FALSE)
  }
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
  apply_overlay(args[[1L]], args[[2L]],
                dirname(normalizePath(script_path, mustWork = TRUE)))
}
