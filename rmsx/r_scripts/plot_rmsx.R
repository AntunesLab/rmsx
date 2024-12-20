# # could add these to the requirements
# library(ggplot2)
# library(viridis)
# library(tidyverse)
# library(readr)
# # library(superheat)
# library(reshape2)
# library(gridExtra)
# library(grid)
# library(cowplot)
#
# # copied from triple_plot_rmsx.R

# Ensure the cluster has the necessary libraries
packages <- c("viridis", "tidyverse")

options(repos = list(CRAN = "https://cloud.r-project.org"))

install_if_not_present <- function(pkg) {
  # Installs a package if it's not already installed, and then loads it.
  #
  # Args:
  #   pkg: A string with the package name to install/load.
  #
  # This function checks if the package is installed. If not, it installs it from CRAN,
  # then loads the library. Useful for ensuring that all dependencies are satisfied
  # before running the main script.
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

sapply(packages, install_if_not_present)

####
# List of required packages
packages <- c("ggplot2", "viridis", "dplyr", "tidyr", "stringr", "readr", "gridExtra", "grid")

sapply(packages, install_if_not_present)


#' @title parse_args
#' @description Parses the command line arguments provided to the script.
#'
#' @return A list with elements:
#'   - csv_path: The path to the RMSX CSV data file.
#'   - rmsd: The path to the RMSD CSV file.
#'   - rmsf: The path to the RMSF CSV file.
#'   - interpolate: Logical flag for raster interpolation.
#'   - triple: Logical flag indicating whether to produce a triple plot (RMSX, RMSD, RMSF).
#'   - palette: The color palette name (for viridis).
#'
#' @details This function reads arguments passed to the script. If no arguments are found,
#' it stops. Otherwise, it prints them for debugging and returns them in a named list.
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    stop("No CSV file path provided!", call. = FALSE)
  }

  print("arg 1"); print(args[1])
  print("arg 2"); print(args[2])
  print("arg 3"); print(args[3])
  print("arg 4"); print(args[4])
  print("arg 5"); print(args[5])
  print("arg 6"); print(args[6])

  list(
    csv_path = args[1],
    rmsd = args[2],
    rmsf = args[3],
    interpolate = ifelse("TRUE" == args[4], TRUE, FALSE),
    triple = ifelse("TRUE" == args[5], TRUE, FALSE),
    palette = args[6]
  )
}

#' @title read_and_summarize_csv
#' @description Reads the RMSX CSV file and prints the count of rows by ChainID.
#'
#' @param csv_path String, path to the RMSX CSV file.
#' @return A data frame (tibble) of the raw RMSX data.
#'
#' @details This function reads the RMSX CSV and groups by ChainID to show how many
#' residues or entries exist per chain. This helps in verifying that the input data
#' is as expected.
read_and_summarize_csv <- function(csv_path) {
  rmsx_raw <- read_csv(csv_path)
  count_by_chainID <- rmsx_raw %>%
    group_by(ChainID) %>%
    summarise(Count = n(), .groups = 'drop')
  print(count_by_chainID)
  rmsx_raw
}

#' @title process_data_by_chain_id
#' @description Processes RMSX data for a specific Chain ID and prepares it for plotting.
#'
#' @param rmsx_raw A data frame (tibble) containing the RMSX data.
#' @param id The chain identifier (e.g., "A").
#' @param csv_path Path to the CSV file containing the RMSX data.
#' @param interpolate Logical, whether to interpolate in geom_raster.
#' @param palette The viridis color palette option to use.
#' @return A ggplot object representing the RMSX plot for the given chain.
#'
#' @details
#' This function extracts data for one specific Chain ID. It determines the simulation length
#' from the filename and calculates the spacing (step_size) between time points.
#'
#' **Key Change and Reasoning:**
#' Instead of labeling time points starting at zero or at the first step, we position the
#' cell centers so that the entire set of cells fits exactly from 0 to sim_len.
#' We do this by starting the time points at `step_size/2` and ending at `sim_len - step_size/2`.
#' This ensures that the left edge of the first cell aligns with time = 0 and the right edge
#' of the last cell aligns with time = sim_len, so no cells are half-cut or misaligned.
process_data_by_chain_id <- function(rmsx_raw, id, csv_path, interpolate, palette) {
  rmsx <- rmsx_raw %>%
    filter(ChainID == id) %>%
    select(-ChainID)

  filename <- basename(csv_path)
  parts <- str_split(filename, "_", simplify = TRUE)
  simulation_time <- as.numeric(parts[ncol(parts) - 1])
  message("Simulation time =", simulation_time)

  sim_len <- simulation_time
  num_time_points <- ncol(rmsx) - 1
  step_size <- sim_len / num_time_points

  # Centers run from step_size/2 to sim_len - step_size/2
  # This makes the first raster cell start at time=0 and the last end exactly at sim_len.
  column_numbers <- seq(step_size/2, sim_len - step_size/2, length.out = num_time_points)

  names(rmsx) <- c("Residue", column_numbers)
  rmsx_long <- pivot_longer(rmsx, cols = -Residue, names_to = "Time_Point", values_to = "RMSF")
  rmsx_long$Time_Point <- as.numeric(rmsx_long$Time_Point)

  plot_rmsx(rmsx_long, interpolate, palette, step_size, sim_len)
}

#' @title plot_rmsx
#' @description Creates a raster plot of RMSX values over time and residue index.
#'
#' @param rmsx_long A long-format data frame with columns: Residue, Time_Point, and RMSF.
#' @param interpolate Logical, whether geom_raster should interpolate the fill values.
#' @param palette The viridis palette option.
#' @param step_size The time increment between frames.
#' @param sim_len The total simulation length in ns.
#' @return A ggplot object of the RMSX raster plot.
#'
#' @details
#' This function plots RMSX values as a heatmap (raster) over time (x-axis) and residue (y-axis).
#'
#' **Key Change and Reasoning:**
#' The `coord_cartesian(xlim = c(0, sim_len))` ensures the plot covers the full simulation time exactly.
#' Because we set the centers of the cells to be from step_size/2 to sim_len - step_size/2, the left
#' edge of the first cell is at 0 and the right edge of the last cell is at sim_len. Thus, the raster
#' perfectly matches the simulation range without extra padding or half-cells at the ends.
plot_rmsx <- function(rmsx_long, interpolate, palette, step_size, sim_len) {
  ggplot(rmsx_long, aes(Time_Point, Residue, fill = RMSF)) +
    geom_raster(interpolate = interpolate) +
    scale_fill_viridis(option = palette) +
    coord_cartesian(xlim = c(0, sim_len)) +
    theme_minimal() +
    theme(legend.position = "left") +
    labs(x = "Time (ns)", y = "Residue (Index)", fill = "RMSX")
}

#' @title save_plot
#' @description Saves the plot to a PNG file.
#'
#' @param rmsx_plot A ggplot object to save.
#' @param csv_path Path to the RMSX CSV file used to generate the plot.
#' @param id The Chain ID, used to name the output file.
#' @return The output file path as a string.
#'
#' @details
#' This function uses ggsave to export the plot to a PNG. The filename is derived from the CSV file name,
#' inserting the chain ID so that multiple chains produce separate files.
save_plot <- function(rmsx_plot, csv_path, id) {
  output_filename <- gsub(pattern = "\\.csv$", replacement = paste("_rmsx_plot_chain_", id, ".png", sep=""), basename(csv_path))
  output_filepath <- file.path(dirname(csv_path), output_filename)
  ggsave(output_filepath, plot = rmsx_plot, width = 10, height = 6, dpi = 200)
  output_filepath
}

#' @title plot_rmsd
#' @description Creates a line plot of RMSD values over frames.
#'
#' @param rmsd A data frame containing RMSD values.
#' @param Frame The frame column from the RMSD data.
#' @param RMSD The RMSD column.
#' @return A ggplot object showing RMSD over frames.
plot_rmsd <- function(rmsd, Frame, RMSD) {
  ggplot(rmsd, aes(x=Frame, y=RMSD)) +
    geom_line() +
    theme_minimal() +
    labs(x = "")
}

#' @title plot_rmsf
#' @description Creates a line plot of RMSF values over residues.
#'
#' @param rmsf_whole_traj A data frame with RMSF values.
#' @param ResidueID The residue ID column.
#' @param RMSF The RMSF column.
#' @return A ggplot object showing RMSF per residue.
plot_rmsf <- function(rmsf_whole_traj, ResidueID, RMSF) {
  ggplot(rmsf_whole_traj, aes(x=ResidueID, y=RMSF)) +
    geom_line() +
    theme_minimal() +
    coord_flip() +
    labs(x = "")
}

#' @title plot_triple
#' @description Arranges RMSD, RMSX, and RMSF plots into a single composite figure.
#'
#' @param rmsx_plot The RMSX heatmap plot.
#' @param rmsd_plot The RMSD line plot.
#' @param rmsf_plot The RMSF line plot.
#' @return A grob containing the arranged plots.
#'
#' @details
#' This function arranges the three plots (RMSD, RMSX, RMSF) into a grid layout for easier comparison.
plot_triple <- function(rmsx_plot, rmsd_plot, rmsf_plot) {
  plot_elements <- list(rmsd_plot, rmsx_plot, rmsf_plot, nullGrob())

  arranged_elements <- arrangeGrob(grobs = plot_elements, layout_matrix=
                                     rbind(c(4,4,1,1,1,1,1,1,4,4,4),
                                           c(4,4,1,1,1,1,1,1,4,4,4),
                                           c(4,2,2,2,2,2,2,2,3,3,3),
                                           c(4,2,2,2,2,2,2,2,3,3,3),
                                           c(4,2,2,2,2,2,2,2,3,3,3),
                                           c(4,2,2,2,2,2,2,2,3,3,3)))
  arranged_elements
}

#' @title main
#' @description The main function executing the analysis and plotting.
#'
#' @details
#' - Parses arguments to get CSV paths and options.
#' - Reads the RMSX data and summarizes it.
#' - For each Chain ID in the data:
#'   1. Processes the data to create RMSX plots with proper time alignment.
#'   2. Saves the RMSX plot.
#'   3. If triple plotting is enabled, also reads RMSD and RMSF data, creates the combined triple plot, and saves it.
#' - Finally, prints out where the plots are saved.
main <- function() {
  args <- parse_args()
  rmsx_raw <- read_and_summarize_csv(args$csv_path)

  for (id in unique(rmsx_raw$ChainID)) {
    rmsx_plot <- process_data_by_chain_id(rmsx_raw, id, args$csv_path, args$interpolate, args$palette)
    save_plot(rmsx_plot, args$csv_path, id)

    if (args$triple == TRUE) {
      rmsd_data <- read_csv(args$rmsd)
      rmsf_data <- read_csv(args$rmsf)

      rmsd_plot <- plot_rmsd(rmsd_data, rmsd_data$Frame, rmsd_data$RMSD)
      rmsf_plot <- plot_rmsf(rmsf_data, rmsf_data$ResidueID, rmsf_data$RMSF)
      triple_plot <- plot_triple(rmsx_plot, rmsd_plot, rmsf_plot)
      save_plot(triple_plot, args$csv_path, id)

      output_file <- "./triple_plot.png"
      png(output_file, width = 800, height = 600)
      grid::grid.draw(triple_plot)
    } else {
      output_file <- "./rmsx_plot.png"
      png(output_file, width = 800, height = 600)
      grid::grid.draw(rmsx_plot)
    }

    dev.off()
    cat("Plot saved to:", normalizePath(output_file), "\n")
    rmsx_plot
  }
}

main()
