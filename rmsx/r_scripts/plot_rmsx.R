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
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

base_needed <- c("ggplot2","viridis","dplyr","tidyr","stringr","readr","gridExtra","grid")
base_missing <- base_needed[!vapply(base_needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(base_missing)) {
  message("RMSX: installing: ", paste(base_missing, collapse = ", "))
  install.packages(base_missing, lib = user_lib,
                   dependencies = c("Depends","Imports","LinkingTo"))
}
suppressPackageStartupMessages(invisible(lapply(base_needed, require, character.only = TRUE)))
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(viridis); library(dplyr); library(tidyr)
  library(stringr); library(readr);   library(gridExtra); library(grid)
})

# ---- argument parsing ------------------------------------------------------------------
# args: 1 csv_path, 2 rmsd_csv, 3 rmsf_csv, 4 interpolate(T/F), 5 triple(T/F),
#       6 palette, 7 manual_min ("" or num), 8 manual_max ("" or num),
#       9 log_transform (T/F), 10 custom_fill_label ("" or text)
#       11 window_check (T/F; optional, default FALSE)         # NEW
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 6) {
    stop("Usage: plot_rmsx.R <rmsx_csv> <rmsd_csv> <rmsf_csv> <interpolate TRUE|FALSE> <triple TRUE|FALSE> <palette> [min] [max] [log TRUE|FALSE] [fill_label] [window_check TRUE|FALSE]",
         call. = FALSE)
  }
  
  manual_min <- if (length(args) >= 7 && nzchar(args[7])) suppressWarnings(as.numeric(args[7])) else NA_real_
  manual_max <- if (length(args) >= 8 && nzchar(args[8])) suppressWarnings(as.numeric(args[8])) else NA_real_
  log_transform <- if (length(args) >= 9 && nzchar(args[9])) identical(toupper(args[9]), "TRUE") else FALSE
  custom_fill_label <- if (length(args) >= 10) args[10] else ""
  window_check <- if (length(args) >= 11 && nzchar(args[11])) identical(toupper(args[11]), "TRUE") else FALSE  # NEW
  
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
    custom_fill_label = custom_fill_label,
    window_check = window_check   # NEW
  )
}

# ---- IO helpers ------------------------------------------------------------------------
read_and_summarize_csv <- function(csv_path) {
  rmsx_raw <- readr::read_csv(csv_path, show_col_types = FALSE)
  by_chain <- rmsx_raw %>% group_by(ChainID) %>% summarise(Count = n(), .groups = "drop")
  print(by_chain)
  rmsx_raw
}

load_mask_metadata <- function(csv_path) {
  mask_path <- file.path(dirname(csv_path), "masked_residues.csv")
  if (!file.exists(mask_path)) {
    return(NULL)
  }

  mask_data <- readr::read_csv(mask_path, show_col_types = FALSE)
  required <- c("ResidueID", "ChainID", "Masked")
  if (!all(required %in% names(mask_data))) {
    stop("Mask metadata is missing required columns: ResidueID, ChainID, Masked", call. = FALSE)
  }

  mask_data %>%
    transmute(
      ResidueID = as.integer(ResidueID),
      ChainID   = as.character(ChainID),
      Masked    = as.logical(Masked)
    )
}

ensure_pattern_support <- function(mask_data) {
  if (is.null(mask_data) || !any(mask_data$Masked, na.rm = TRUE)) {
    return(FALSE)
  }

  if (getRversion() < "4.1.0") {
    stop(
      "Masked heatmap rendering requires R 4.1+ because RMSX uses ggpattern for hatch overlays.",
      call. = FALSE
    )
  }

  pattern_needed <- c("ggpattern")
  pattern_missing <- pattern_needed[!vapply(pattern_needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(pattern_missing)) {
    message("RMSX: installing masked-plot support: ", paste(pattern_missing, collapse = ", "))
    install.packages(pattern_missing, lib = user_lib,
                     dependencies = c("Depends","Imports","LinkingTo"))
  }
  suppressPackageStartupMessages(invisible(lapply(pattern_needed, require, character.only = TRUE)))
  TRUE
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
  
  base_plot <- ggplot(rmsx_long, aes(Time_Point, Residue, fill = RMSF)) +
    geom_raster(interpolate = interpolate) +
    fill_scale +
    coord_cartesian(xlim = c(0, sim_len)) +
    theme_minimal() +
    theme(legend.position = "left") +
    labs(x = "Time (ns)", y = "Residue (Index)", fill = fill_label)

  if ("Masked" %in% names(rmsx_long) && any(rmsx_long$Masked, na.rm = TRUE)) {
    masked_tiles <- rmsx_long %>% filter(Masked)
    base_plot <- base_plot +
      ggpattern::geom_tile_pattern(
        data = masked_tiles,
        inherit.aes = FALSE,
        aes(x = Time_Point, y = Residue, pattern = "Masked residues"),
        fill = NA,
        color = NA,
        width = step_size,
        height = 1,
        pattern_fill = "black",
        pattern_colour = "black",
        pattern_alpha = 0.55,
        pattern_angle = 45,
        pattern_density = 0.30,
        pattern_spacing = 0.02,
        pattern_key_scale_factor = 0.6
      ) +
      ggpattern::scale_pattern_manual(values = c("Masked residues" = "stripe"), name = NULL)
  }

  base_plot
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

# ---- NEW: per-slice mean RMSX plot -----------------------------------------------------
# Given the long-format RMSX (one row per residue × time),
# compute the mean across residues for each time point and plot it as a line.
plot_rmsx_mean <- function(rmsx_long) {
  mean_df <- rmsx_long %>%
    group_by(Time_Point) %>%
    summarise(Mean_RMSX = mean(RMSF, na.rm = TRUE), .groups = "drop")
  
  ggplot(mean_df, aes(x = Time_Point, y = Mean_RMSX)) +
    geom_line() +
    theme_minimal() +
    labs(x = "", y = "Mean RMSX")
}

# ---- NEW: compute slice-mean RMSD aligned to RMSX slices ------------

compute_mean_rmsd_per_slice <- function(rmsd, num_slices, sim_len) {
  if (num_slices <= 0 || sim_len <= 0) {
    stop("compute_mean_rmsd_per_slice(): num_slices and sim_len must be positive.")
  }
  
  # Use "Time" in ps if available, otherwise spread frames evenly over [0, sim_len]
  if ("Time" %in% names(rmsd) && all(is.finite(rmsd$Time))) {
    t_ns <- rmsd$Time / 1000.0
  } else {
    n <- nrow(rmsd)
    if (n <= 1) {
      t_ns <- rep(0, n)
    } else {
      t_ns <- seq(0, sim_len, length.out = n)
    }
  }
  
  df <- rmsd %>%
    mutate(Time_ns = t_ns)
  
  step <- sim_len / num_slices
  
  df <- df %>%
    mutate(
      Slice = floor(Time_ns / step) + 1L,
      Slice = pmin(Slice, num_slices)
    )
  
  slice_df <- df %>%
    group_by(Slice) %>%
    summarise(
      Slice_center = (Slice - 0.5) * step,
      Mean_RMSD    = mean(RMSD, na.rm = TRUE),
      .groups = "drop"
    )
  
  slice_df
}

# ---- NEW: slice-mean RMSD plot, aligned to RMSX slices --------------
# Given the RMSD time series and the RMSX time grid (sim_len, num_slices),
# compute mean RMSD per slice and plot it vs slice-center time (ns).
plot_rmsd_sliced <- function(rmsd, num_slices, sim_len) {
  slice_df <- compute_mean_rmsd_per_slice(rmsd, num_slices, sim_len)
  
  ggplot(slice_df, aes(x = Slice_center, y = Mean_RMSD)) +
    geom_line() +                         # line only, no points
    theme_minimal() +
    coord_cartesian(xlim = c(0, sim_len)) +
    labs(x = "", y = "Mean RMSD")
}

# ---- NEW: stacked RMSD + slice-mean RMSD + mean RMSX + heatmap -----------
plot_stacked_rmsd_mean_rmsx_heatmap <- function(rmsd_plot,
                                                rmsd_sliced_plot,
                                                rmsx_mean_plot,
                                                rmsx_plot) {
  # 1 = RMSD (top)
  # 2 = slice-mean RMSD
  # 3 = mean RMSX
  # 4 = RMSX heatmap (bottom)
  # 5 = whitespace / nullGrob()
  arrangeGrob(
    grobs = list(
      rmsd_plot,        # 1
      rmsd_sliced_plot, # 2
      rmsx_mean_plot,   # 3
      rmsx_plot,        # 4
      nullGrob()        # 5 (whitespace)
    ),
    layout_matrix = rbind(
      c(5,5,1,1,1,1,1,1,5,5,5),
      c(5,5,1,1,1,1,1,1,5,5,5),
      c(5,5,1,1,1,1,1,1,5,5,5),
      c(5,5,1,1,1,1,1,1,5,5,5),
      c(5,5,1,1,1,1,1,1,5,5,5),
      c(5,5,1,1,1,1,1,1,5,5,5),
      c(5,5,2,2,2,2,2,2,5,5,5),
      c(5,5,2,2,2,2,2,2,5,5,5),
      c(5,5,2,2,2,2,2,2,5,5,5),
      c(5,5,2,2,2,2,2,2,5,5,5),
      c(5,5,3,3,3,3,3,3,5,5,5),
      c(5,5,3,3,3,3,3,3,5,5,5),
      c(5,5,3,3,3,3,3,3,5,5,5),
      c(5,5,3,3,3,3,3,3,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5),
      c(5,4,4,4,4,4,4,4,5,5,5)
    )
  )
}

save_plot <- function(plot_obj, csv_path, id) {
  outfile <- sub("\\.csv$", paste0("_rmsx_plot_chain_", id, ".png"), basename(csv_path))
  outpath <- file.path(dirname(csv_path), outfile)
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
  
  filename <- basename(csv_path)
  parts <- stringr::str_split(filename, "_", simplify = TRUE)
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
  rmsx_long <- tidyr::pivot_longer(rmsx, cols = -Residue, names_to = "Time_Point", values_to = "RMSF")
  rmsx_long$Time_Point <- as.numeric(rmsx_long$Time_Point)

  mask_data <- load_mask_metadata(csv_path)
  if (!is.null(mask_data)) {
    mask_chain <- mask_data %>%
      filter(ChainID == as.character(id)) %>%
      transmute(Residue = ResidueID, Masked = Masked)
    if (nrow(mask_chain)) {
      rmsx_long <- rmsx_long %>%
        left_join(mask_chain, by = "Residue") %>%
        mutate(Masked = if_else(is.na(Masked), FALSE, Masked))
    } else {
      rmsx_long$Masked <- FALSE
    }
  } else {
    rmsx_long$Masked <- FALSE
  }
  
  # Return both the heatmap, the long data, and the time grid info
  list(
    heatmap         = plot_rmsx(rmsx_long, interpolate, palette, step_size, sim_len,
                                manual_min, manual_max, log_transform, custom_fill_label),
    long            = rmsx_long,
    sim_len         = sim_len,
    num_time_points = num_time_points,
    step_size       = step_size
  )
}

# ---- NEW: compute correlations for window_check ---------------------
# Uses the wide-format rmsx_raw (with ResidueID, ChainID, slice_* columns),
# plus RMSD and RMSF CSVs, to compute:
#   1) Residue-level correlation: RMSF vs mean RMSX (across slices)
#   2) Global correlation: mean RMSX(slice) vs mean RMSD(slice)
compute_window_check_correlations <- function(rmsx_raw, chain_id, rmsd_data, rmsf_data) {
  if (!("ResidueID" %in% names(rmsx_raw))) {
    warning("RMSX CSV missing 'ResidueID'; skipping correlation calculation.")
    return(invisible(NULL))
  }
  if (!("ChainID" %in% names(rmsx_raw))) {
    warning("RMSX CSV missing 'ChainID'; skipping correlation calculation.")
    return(invisible(NULL))
  }

  rmsx_chain <- rmsx_raw %>% filter(ChainID == chain_id)
  if (nrow(rmsx_chain) == 0) {
    warning("No rows for ChainID ", chain_id, " in RMSX CSV; skipping correlations.")
    return(invisible(NULL))
  }

  slice_cols <- setdiff(names(rmsx_chain), c("ResidueID", "ChainID"))
  if (length(slice_cols) < 1L) {
    warning("No slice columns found for ChainID ", chain_id, "; skipping correlations.")
    return(invisible(NULL))
  }

  # 1) Residue-level: RMSF vs mean RMSX ------------------------------
  mean_rmsx_res <- rmsx_chain %>%
    mutate(mean_RMSX = rowMeans(across(all_of(slice_cols)), na.rm = TRUE)) %>%
    select(ResidueID, mean_RMSX)

  rmsf_chain <- rmsf_data
  if ("ChainID" %in% names(rmsf_chain) && !is.na(chain_id)) {
    rmsf_chain <- rmsf_chain %>% filter(ChainID == chain_id)
  }

  r_res  <- NA_real_
  r2_res <- NA_real_
  if (all(c("ResidueID", "RMSF") %in% names(rmsf_chain))) {
    res_join <- rmsf_chain %>%
      select(ResidueID, RMSF) %>%
      inner_join(mean_rmsx_res, by = "ResidueID")

    if (nrow(res_join) > 0) {
      r_res  <- suppressWarnings(cor(res_join$RMSF, res_join$mean_RMSX, use = "complete.obs"))
      r2_res <- r_res^2
    } else {
      warning("No overlapping ResidueID between RMSF and RMSX for ChainID ", chain_id)
    }
  } else {
    warning("RMSF CSV missing 'ResidueID' or 'RMSF' for residue-level correlation.")
  }

  # 2) Global-level: mean RMSX per slice vs mean RMSD per slice -------
  means <- colMeans(rmsx_chain[, slice_cols, drop = FALSE], na.rm = TRUE)
  mean_rmsx_slice <- data.frame(
    Slice     = seq_along(means),
    mean_RMSX = as.numeric(means),
    stringsAsFactors = FALSE
  )
  n_slices <- nrow(mean_rmsx_slice)

  r_glob  <- NA_real_
  r2_glob <- NA_real_
  if (!("RMSD" %in% names(rmsd_data))) {
    warning("RMSD CSV missing 'RMSD'; skipping global correlation.")
  } else {
    # Build a time axis from 'Time' (if present) or 'Frame'
    if ("Time" %in% names(rmsd_data) && all(is.finite(rmsd_data$Time))) {
      t <- rmsd_data$Time
    } else if ("Frame" %in% names(rmsd_data) && all(is.finite(rmsd_data$Frame))) {
      t <- rmsd_data$Frame
    } else {
      warning("RMSD CSV must have either 'Time' or 'Frame' for global correlation.")
      t <- NULL
    }

    if (!is.null(t)) {
      t <- t - min(t, na.rm = TRUE)
      t_max <- max(t, na.rm = TRUE)
      if (!is.finite(t_max) || t_max == 0) t_max <- 1

      slice_len <- t_max / n_slices
      slice_index <- floor(t / slice_len) + 1L
      slice_index[slice_index > n_slices] <- n_slices

      df <- data.frame(
        Slice = slice_index,
        RMSD  = rmsd_data$RMSD
      )

      mean_rmsd_slice <- df %>%
        group_by(Slice) %>%
        summarise(mean_RMSD = mean(RMSD, na.rm = TRUE),
                  .groups = "drop")

      all_slices <- data.frame(Slice = seq_len(n_slices))
      mean_rmsd_slice <- all_slices %>%
        left_join(mean_rmsd_slice, by = "Slice")

      slice_join <- mean_rmsx_slice %>%
        inner_join(mean_rmsd_slice, by = "Slice")

      if (nrow(slice_join) > 0) {
        r_glob  <- suppressWarnings(cor(slice_join$mean_RMSX, slice_join$mean_RMSD,
                                        use = "complete.obs"))
        r2_glob <- r_glob^2
      } else {
        warning("No overlapping slices when joining mean RMSX and mean RMSD for ChainID ",
                chain_id)
      }
    }
  }
}

# ---- NEW: compute correlations for window_check ---------------------
# Uses the wide-format rmsx_raw (with ResidueID, ChainID, slice_* columns),
# plus RMSD and RMSF CSVs, to compute:
#   1) Residue-level correlation: RMSF vs mean RMSX (across slices)
#   2) Global correlation: per-frame RMSD vs mean RMSX of that frame's slice
# compute_window_check_correlations <- function(rmsx_raw, chain_id, rmsd_data, rmsf_data) {
#   if (!("ResidueID" %in% names(rmsx_raw))) {
#     warning("RMSX CSV missing 'ResidueID'; skipping correlation calculation.")
#     return(invisible(NULL))
#   }
#   if (!("ChainID" %in% names(rmsx_raw))) {
#     warning("RMSX CSV missing 'ChainID'; skipping correlation calculation.")
#     return(invisible(NULL))
#   }
#   
#   rmsx_chain <- rmsx_raw %>% filter(ChainID == chain_id)
#   if (nrow(rmsx_chain) == 0) {
#     warning("No rows for ChainID ", chain_id, " in RMSX CSV; skipping correlations.")
#     return(invisible(NULL))
#   }
#   
#   slice_cols <- setdiff(names(rmsx_chain), c("ResidueID", "ChainID"))
#   if (length(slice_cols) < 1L) {
#     warning("No slice columns found for ChainID ", chain_id, "; skipping correlations.")
#     return(invisible(NULL))
#   }
#   
#   # 1) Residue-level: RMSF vs mean RMSX ------------------------------
#   mean_rmsx_res <- rmsx_chain %>%
#     mutate(mean_RMSX = rowMeans(across(all_of(slice_cols)), na.rm = TRUE)) %>%
#     select(ResidueID, mean_RMSX)
#   
#   rmsf_chain <- rmsf_data
#   if ("ChainID" %in% names(rmsf_chain) && !is.na(chain_id)) {
#     rmsf_chain <- rmsf_chain %>% filter(ChainID == chain_id)
#   }
#   
#   r_res  <- NA_real_
#   r2_res <- NA_real_
#   if (all(c("ResidueID", "RMSF") %in% names(rmsf_chain))) {
#     res_join <- rmsf_chain %>%
#       select(ResidueID, RMSF) %>%
#       inner_join(mean_rmsx_res, by = "ResidueID")
#     
#     if (nrow(res_join) > 0) {
#       r_res  <- suppressWarnings(cor(res_join$RMSF, res_join$mean_RMSX, use = "complete.obs"))
#       r2_res <- r_res^2
#     } else {
#       warning("No overlapping ResidueID between RMSF and RMSX for ChainID ", chain_id)
#     }
#   } else {
#     warning("RMSF CSV missing 'ResidueID' or 'RMSF' for residue-level correlation.")
#   }
#   
#   # 2) Global-level: per-frame RMSD vs mean RMSX of that frame's slice ----
#   means <- colMeans(rmsx_chain[, slice_cols, drop = FALSE], na.rm = TRUE)
#   mean_rmsx_slice <- data.frame(
#     Slice     = seq_along(means),
#     mean_RMSX = as.numeric(means),
#     stringsAsFactors = FALSE
#   )
#   n_slices <- nrow(mean_rmsx_slice)
#   
#   r_glob  <- NA_real_
#   r2_glob <- NA_real_
#   if (!("RMSD" %in% names(rmsd_data))) {
#     warning("RMSD CSV missing 'RMSD'; skipping global correlation.")
#   } else {
#     # Build a time axis from 'Time' (if present) or 'Frame'
#     if ("Time" %in% names(rmsd_data) && all(is.finite(rmsd_data$Time))) {
#       t <- rmsd_data$Time
#     } else if ("Frame" %in% names(rmsd_data) && all(is.finite(rmsd_data$Frame))) {
#       t <- rmsd_data$Frame
#     } else {
#       warning("RMSD CSV must have either 'Time' or 'Frame' for global correlation.")
#       t <- NULL
#     }
#     
#     if (!is.null(t)) {
#       t <- t - min(t, na.rm = TRUE)
#       t_max <- max(t, na.rm = TRUE)
#       if (!is.finite(t_max) || t_max == 0) t_max <- 1
#       
#       slice_len   <- t_max / n_slices
#       slice_index <- floor(t / slice_len) + 1L
#       slice_index[slice_index > n_slices] <- n_slices
#       
#       # For each frame, attach the mean RMSX of its slice
#       df <- data.frame(
#         Slice = slice_index,
#         RMSD  = rmsd_data$RMSD
#       ) %>%
#         left_join(mean_rmsx_slice, by = "Slice")
#       
#       valid <- is.finite(df$RMSD) & is.finite(df$mean_RMSX)
#       if (any(valid)) {
#         r_glob  <- suppressWarnings(cor(df$RMSD[valid], df$mean_RMSX[valid],
#                                         use = "complete.obs"))
#         r2_glob <- r_glob^2
#       } else {
#         warning("No finite pairs for global correlation (RMSD vs mean RMSX).")
#       }
#     }
#   }
#   
#   # --- print summary to console (CLI or notebook) --------------------
#   cat("------------------------------------------------------------\n")
#   cat("Window-check correlations for ChainID ", chain_id, ":\n", sep = "")
#   if (!is.na(r_res)) {
#     cat(sprintf("  Residue-level (RMSF vs mean RMSX):            r = %.3f, R^2 = %.3f\n",
#                 r_res, r2_res))
#   } else {
#     cat("  Residue-level (RMSF vs mean RMSX):            not available.\n")
#   }
#   if (!is.na(r_glob)) {
#     cat(sprintf("  Global-level (RMSD vs mean RMSX(slice)):     r = %.3f, R^2 = %.3f\n",
#                 r_glob, r2_glob))
#   } else {
#     cat("  Global-level (RMSD vs mean RMSX(slice)):     not available.\n")
#   }
#   cat("------------------------------------------------------------\n")
#   
#   invisible(list(
#     r_residue  = r_res,
#     r2_residue = r2_res,
#     r_global   = r_glob,
#     r2_global  = r2_glob,
#     n_slices   = n_slices
#   ))
# }
# 
#   
#   # --- print summary to console (CLI or notebook) --------------------
#   cat("------------------------------------------------------------\n")
#   cat("Window-check correlations for ChainID ", chain_id, ":\n", sep = "")
#   if (!is.na(r_res)) {
#     cat(sprintf("  Residue-level (RMSF vs mean RMSX):      r = %.3f, R^2 = %.3f\n",
#                 r_res, r2_res))
#   } else {
#     cat("  Residue-level (RMSF vs mean RMSX):      not available.\n")
#   }
#   if (!is.na(r_glob)) {
#     cat(sprintf("  Global-level (mean RMSX vs mean RMSD):  r = %.3f, R^2 = %.3f\n",
#                 r_glob, r2_glob))
#   } else {
#     cat("  Global-level (mean RMSX vs mean RMSD):  not available.\n")
#   }
#   cat("------------------------------------------------------------\n")
#   
#   invisible(list(
#     r_residue  = r_res,
#     r2_residue = r2_res,
#     r_global   = r_glob,
#     r2_global  = r2_glob,
#     n_slices   = n_slices
#   ))
# }

# ---- main ------------------------------------------------------------------------------
main <- function() {
  args <- parse_args()
  rmsx_raw <- read_and_summarize_csv(args$csv_path)
  mask_data <- load_mask_metadata(args$csv_path)
  ensure_pattern_support(mask_data)
  
  if (isTRUE(args$log_transform)) {
    message("Log transform requested: assuming CSV already contains log-scaled values.")
  }
  
  for (id in unique(rmsx_raw$ChainID)) {
    res <- process_data_by_chain_id(
      rmsx_raw, id, args$csv_path,
      args$interpolate, args$palette,
      args$manual_min, args$manual_max,
      args$log_transform, args$custom_fill_label
    )
    rmsx_plot <- res$heatmap
    rmsx_long <- res$long
    
    if (isTRUE(args$triple)) {
      rmsd_data <- readr::read_csv(args$rmsd, show_col_types = FALSE)
      rmsf_data <- readr::read_csv(args$rmsf, show_col_types = FALSE)
      
      if (isTRUE(args$window_check)) {
        # NEW: stacked RMSD + slice-mean RMSD + mean RMSX + heatmap (no RMSF)
        rmsd_plot        <- plot_rmsd(rmsd_data)
        rmsd_sliced_plot <- plot_rmsd_sliced(
          rmsd        = rmsd_data,
          num_slices  = res$num_time_points,
          sim_len     = res$sim_len
        )
        rmsx_mean_plot   <- plot_rmsx_mean(rmsx_long)
        
        stacked_plot <- plot_stacked_rmsd_mean_rmsx_heatmap(
          rmsd_plot,
          rmsd_sliced_plot,
          rmsx_mean_plot,
          rmsx_plot
        )
        save_plot(stacked_plot, args$csv_path, id)
        
        # NEW: compute and print correlations as part of window_check
        compute_window_check_correlations(
          rmsx_raw   = rmsx_raw,
          chain_id   = id,
          rmsd_data  = rmsd_data,
          rmsf_data  = rmsf_data
        )
      } else {
        # Original triple layout (RMSD + RMSX heatmap + RMSF)
        triple_plot <- plot_triple(rmsx_plot, plot_rmsd(rmsd_data), plot_rmsf(rmsf_data))
        save_plot(triple_plot, args$csv_path, id)   # overwrites heatmap filename
      }
    } else {
      save_plot(rmsx_plot, args$csv_path, id)
    }
  }
}

invisible(main())
