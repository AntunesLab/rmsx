#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------------------
# RMSX heatmaps (and optional triple plots with RMSD/RMSF)
# - Minimal cross-platform package bootstrap (no sudo; installs to user library)
# - Non-interactive: sets CRAN repo and suppresses prompts
# - Same CLI interface as before
# - Saves one PNG per chain: "<csv_basename>_rmsx_plot_chain_<ID>.png"
#   (If triple=TRUE, the triple composite overwrites that filename.)
# ----------------------------------------------------------------------------------------


# ---- minimal package bootstrap ---------------------------------------------------------
if (Sys.getenv("RMSX_NO_AUTO_INSTALL", "0") != "1") {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (user_lib == "") user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))

  # Prefer a repo with binaries when possible; CRAN as fallback
  if (.Platform$OS.type == "unix") {
    options(repos = c(RSPM = "https://packagemanager.posit.co/cran/latest",
                      CRAN = "https://cloud.r-project.org"))
  } else {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }

  needed <- c("ggplot2","viridis","dplyr","tidyr","stringr","readr","gridExtra","grid")
  miss <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(miss)) {
    message("RMSX: installing missing R packages into: ", normalizePath(user_lib, winslash = "/"))
    message("RMSX: missing packages: ", paste(miss, collapse = ", "))

    # Speed up source builds
    if (is.na(getOption("Ncpus"))) {
      nc <- tryCatch(parallel::detectCores(), error = function(e) 1)
      options(Ncpus = max(1, as.integer(nc)))
    }

    t0 <- proc.time()[["elapsed"]]
    ok <- TRUE
    tryCatch(
      install.packages(miss, lib = user_lib, dependencies = TRUE, quiet = FALSE),
      error = function(e) {
        ok <<- FALSE
        message("RMSX: package auto-install failed (no internet or blocked repo?).")
        message("      Preinstall the above packages or set RMSX_NO_AUTO_INSTALL=1 to skip.")
        message("      Details: ", conditionMessage(e))
      }
    )
    t1 <- proc.time()[["elapsed"]]
    if (ok) message(sprintf("RMSX: package install finished in %.1f sec.", t1 - t0))
  } else {
    message("RMSX: all required R packages are already installed.")
  }
}

suppressPackageStartupMessages({
  library(ggplot2); library(viridis); library(dplyr); library(tidyr)
  library(stringr); library(readr);   library(gridExtra); library(grid)
})



# ---- argument parsing ------------------------------------------------------------------
# args: 1 csv_path, 2 rmsd_csv, 3 rmsf_csv, 4 interpolate(T/F), 5 triple(T/F),
#       6 palette, 7 manual_min ("" or num), 8 manual_max ("" or num),
#       9 log_transform (T/F), 10 custom_fill_label ("" or text)
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 6) {
    stop("Usage: plot_rmsx.R <rmsx_csv> <rmsd_csv> <rmsf_csv> <interpolate TRUE|FALSE> <triple TRUE|FALSE> <palette> [min] [max] [log TRUE|FALSE] [fill_label]",
         call. = FALSE)
  }

  # Optional numeric min/max
  manual_min <- if (length(args) >= 7 && nzchar(args[7])) suppressWarnings(as.numeric(args[7])) else NA_real_
  manual_max <- if (length(args) >= 8 && nzchar(args[8])) suppressWarnings(as.numeric(args[8])) else NA_real_

  # Optional log flag
  log_transform <- if (length(args) >= 9 && nzchar(args[9])) identical(toupper(args[9]), "TRUE") else FALSE

  # Optional custom fill label
  custom_fill_label <- if (length(args) >= 10) args[10] else ""

  list(
    csv_path   = args[1],
    rmsd       = args[2],
    rmsf       = args[3],
    interpolate= identical(toupper(args[4]), "TRUE"),
    triple     = identical(toupper(args[5]), "TRUE"),
    palette    = args[6],
    manual_min = manual_min,
    manual_max = manual_max,
    log_transform = log_transform,
    custom_fill_label = custom_fill_label
  )
}

# ---- IO helpers ------------------------------------------------------------------------
read_and_summarize_csv <- function(csv_path) {
  rmsx_raw <- readr::read_csv(csv_path, show_col_types = FALSE)
  by_chain <- rmsx_raw %>% group_by(ChainID) %>% summarise(Count = n(), .groups = "drop")
  print(by_chain)
  rmsx_raw
}

# ---- plotting helpers ------------------------------------------------------------------
plot_rmsx <- function(rmsx_long, interpolate, palette, step_size, sim_len,
                      manual_min, manual_max, log_transform = FALSE, custom_fill_label = "") {

  fill_scale <- if (!is.na(manual_min) && !is.na(manual_max)) {
    scale_fill_viridis(option = palette, limits = c(manual_min, manual_max))
  } else {
    scale_fill_viridis(option = palette)
  }

  fill_label <- if (nzchar(custom_fill_label)) {
    custom_fill_label
  } else if (isTRUE(log_transform)) {
    "Log-\nScaled\nRMSX"
  } else {
    "RMSX (Å)"
  }

  ggplot(rmsx_long, aes(Time_Point, Residue, fill = RMSF)) +
    geom_raster(interpolate = interpolate) +
    fill_scale +
    coord_cartesian(xlim = c(0, sim_len)) +
    theme_minimal() +
    theme(legend.position = "left") +
    labs(x = "Time (ns)", y = "Residue (Index)", fill = fill_label)
}

plot_rmsd <- function(rmsd) {
  ggplot(rmsd, aes(x = Frame, y = RMSD)) +
    geom_line() +
    theme_minimal() +
    labs(x = "", y = "RMSD (Å)", title = "")
}

plot_rmsf <- function(rmsf_whole_traj) {
  ggplot(rmsf_whole_traj, aes(x = ResidueID, y = RMSF)) +
    geom_line() +
    theme_minimal() +
    coord_flip() +
    labs(x = "", y = "RMSF (Å)")
}

plot_triple <- function(rmsx_plot, rmsd_plot, rmsf_plot) {
  arrangeGrob(
    grobs = list(rmsd_plot, rmsx_plot, rmsf_plot, nullGrob()),
    layout_matrix = rbind(
      c(4,4,1,1,1,1,1,1,4,4,4),
      c(4,4,1,1,1,1,1,1,4,4,4),
      c(4,2,2,2,2,2,2,2,3,3,3),
      c(4,2,2,2,2,2,2,2,3,3,3),
      c(4,2,2,2,2,2,2,2,3,3,3),
      c(4,2,2,2,2,2,2,2,3,3,3)
    )
  )
}

save_plot <- function(plot_obj, csv_path, id) {
  outfile <- sub("\\.csv$", paste0("_rmsx_plot_chain_", id, ".png"), basename(csv_path))
  outpath <- file.path(dirname(csv_path), outfile)
  # ggsave can handle both ggplot objects and grobs
  ggsave(outpath, plot = plot_obj, width = 10, height = 6, dpi = 200)
  message("Saved: ", normalizePath(outpath, winslash = "/"))
  outpath
}

# ---- per-chain processing --------------------------------------------------------------
process_data_by_chain_id <- function(rmsx_raw, id, csv_path, interpolate, palette,
                                     manual_min, manual_max, log_transform, custom_fill_label) {
  rmsx <- rmsx_raw %>%
    filter(ChainID == id) %>%
    select(-ChainID)

  # Filename pattern: "..._<length_ns>_ns.csv" (e.g., rmsx_mon_sys_0.015_ns.csv)
  filename <- basename(csv_path)
  parts <- stringr::str_split(filename, "_", simplify = TRUE)
  # The numeric length is the penultimate token before "ns.csv"
  # Example tokens: [rmsx][mon][sys][0.015][ns.csv] -> take index ncol(parts)-1
  simulation_time <- suppressWarnings(as.numeric(parts[ncol(parts) - 1]))
  if (!is.finite(simulation_time)) {
    stop("Could not parse simulation time from filename: ", filename)
  }
  message("Simulation time = ", simulation_time)

  sim_len <- simulation_time
  num_time_points <- ncol(rmsx) - 1
  if (num_time_points <= 0) stop("No time-slice columns found in CSV for ChainID ", id)

  step_size <- sim_len / num_time_points
  centers <- seq(step_size / 2, sim_len - step_size / 2, length.out = num_time_points)

  names(rmsx) <- c("Residue", centers)
  rmsx_long <- pivot_longer(rmsx, cols = -Residue, names_to = "Time_Point", values_to = "RMSF")
  rmsx_long$Time_Point <- as.numeric(rmsx_long$Time_Point)

  plot_rmsx(rmsx_long, interpolate, palette, step_size, sim_len,
            manual_min, manual_max, log_transform, custom_fill_label)
}

# ---- main ------------------------------------------------------------------------------
main <- function() {
  args <- parse_args()
  rmsx_raw <- read_and_summarize_csv(args$csv_path)

  if (isTRUE(args$log_transform)) {
    message("Log transform requested: assuming CSV already contains log-scaled values.")
  }

  for (id in unique(rmsx_raw$ChainID)) {
    rmsx_plot <- process_data_by_chain_id(
      rmsx_raw, id, args$csv_path,
      args$interpolate, args$palette,
      args$manual_min, args$manual_max,
      args$log_transform, args$custom_fill_label
    )

    # If triple, build and save the composite; otherwise save heatmap
    if (isTRUE(args$triple)) {
      rmsd_data <- readr::read_csv(args$rmsd, show_col_types = FALSE)
      rmsf_data <- readr::read_csv(args$rmsf, show_col_types = FALSE)
      triple_plot <- plot_triple(rmsx_plot, plot_rmsd(rmsd_data), plot_rmsf(rmsf_data))
      save_plot(triple_plot, args$csv_path, id)   # overwrites heatmap filename
    } else {
      save_plot(rmsx_plot, args$csv_path, id)
    }
  }
}

invisible(main())

