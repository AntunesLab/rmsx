# Plotting dependencies used by RMSX, including the hatch overlay for masked
# residues.  This script is intended for an explicit first-run setup step.
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

pkgs <- c("ggplot2", "viridis", "dplyr", "tidyr", "stringr", "readr", "gridExtra", "ggpattern")
need <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(need)) {
  message("Installing RMSX plotting support (first run only): ", paste(need, collapse = ", "))
  install.packages(
    need,
    lib = user_lib,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

cat("RMSX plotting support is ready.\n")
