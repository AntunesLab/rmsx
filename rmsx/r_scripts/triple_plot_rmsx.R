# # could add thees to the requirements
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
# 
# # Ensure the cluster has the necessary libraries
# packages <- c("viridis", "tidyverse")
# 
# options(repos = list(CRAN = "https://cloud.r-project.org"))
# 
# install_if_not_present <- function(pkg) {
#   if (!require(pkg, character.only = TRUE)) {
#     install.packages(pkg, dependencies = TRUE)
#     library(pkg, character.only = TRUE)
#   }
# }
# 
# sapply(packages, install_if_not_present)
# 
# # Function to parse command line arguments

####

# List of required packages
packages <- c("ggplot2", "viridis", "dplyr", "tidyr", "stringr", "readr", "gridExtra", "grid")

# Function to install packages if not already installed
install_if_not_present <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Install required packages
sapply(packages, install_if_not_present)


  
  # Simulating command line arguments
  # args <- character(3)
  # args[1] <- "/Users/finn/Desktop/RMSX_Demo_files_mac/ubq_ww_pcv_example/rmsx_ubq_ww_pcv_0.00196_ns.csv"
  # args[2] <- "TRUE"
  # args[3] <- "/Users/finn/Desktop/RMSX_Demo_files_mac/ubq_ww_pcv_example/rmsd.csv"
  # args[4] <- "/Users/finn/Desktop/RMSX_Demo_files_mac/ubq_ww_pcv_example/rmsf.csv"
  #########################################

# # Function to parse command line arguments

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    stop("No CSV file path provided!", call. = FALSE)
  }
  # changed it to 1 or less to make interpolate optional, 
  # list(csv_path = args[1], interpolate = ifelse("TRUE" == args[2], FALSE), rmsd = args[3], rmsf =args[4])
  print("arg 1")
  print(args[1])
  print("arg 2")
  print(args[2])
  print("arg 3")
  print(args[3])
  print("arg 4")
  print(args[4])
  print("arg 5")
  print(args[5])

  list(csv_path = args[1], 
       rmsd = args[2], 
       rmsf = args[3],
       interpolate = ifelse("TRUE" == args[4], TRUE, FALSE), 
       triple = ifelse("TRUE" == args[5], TRUE, FALSE) )
  }


# Function to read and summarize the CSV data
read_and_summarize_csv <- function(csv_path) {
  rmsx_raw <- read_csv(csv_path)
  count_by_chainID <- rmsx_raw %>%
    group_by(ChainID) %>%
    summarise(Count = n(), .groups = 'drop')
  
  print(count_by_chainID)
  rmsx_raw
}

# Function to process data by Chain ID
process_data_by_chain_id <- function(rmsx_raw, id, csv_path, interpolate) {
  rmsx <- rmsx_raw %>%
    filter(ChainID == id) %>%
    select(-ChainID)
  
  filename <- basename(csv_path)
  parts <- str_split(filename, "_", simplify = TRUE)
  simulation_time <- as.numeric(parts[ncol(parts) - 1])
  print(paste("Simulation time =", simulation_time))
  
  sim_len = simulation_time
  step_size <- sim_len / (ncol(rmsx) - 1)
  column_numbers <- seq(from = step_size, to = sim_len, by = step_size)
  names(rmsx) <- c("Residue", column_numbers)
  
  rmsx_long <- pivot_longer(rmsx, cols = -Residue, names_to = "Time_Point", values_to = "RMSF")
  rmsx_long$Time_Point <- as.numeric(rmsx_long$Time_Point)
  
  plot_rmsx(rmsx_long, interpolate)
}

# Function to plot data
plot_rmsx <- function(rmsx_long, interpolate) {
  
  rmsx_plot <- ggplot(rmsx_long, aes(Time_Point, Residue, fill=RMSF)) +
        # geom_point(aes(x=0, y=0), colour="blue") +

    geom_raster(interpolate = interpolate) +
    # facet_wrap(~Chain) +
    scale_fill_viridis(option = "plasma") +
    # it is suprisingly hard to get the y axis aligned at 0. this is a bit of a hack but it works:
    coord_cartesian(xlim=c(0-(rmsx_long$Time_Point[2] - rmsx_long$Time_Point[1])/2,max(rmsx_long$Time_Point))) + # idea  add a space equal to 0 - distance from start to first portion this helps align it to the RMSD plot above, basically just subtract rmsx_long$Time_Point[2] - rmsx_long$Time_Point[1] to get time diff, then take half of it to align it.
  
    # scale_fill_viridis() +
    # theme(aspect.ratio = .1) + # adjust the aspect ratio as needed
    theme_minimal() +
    theme(legend.position = "left")  + # Move legend to the left

    # theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points")) +  # Ensure margins match
    labs(x = "Time", y = "Residue", fill= "RMSX")
  rmsx_plot
  # 
  
  
  rmsx_plot
}

# Function to save the plot
save_plot <- function(rmsx_plot, csv_path, id) {
  output_filename <- gsub(pattern = "\\.csv$", replacement = paste("_rmsx_plot_chain_", id, ".png", sep=""), basename(csv_path))
  output_filepath <- file.path(dirname(csv_path), output_filename)
  ggsave(output_filepath, plot = rmsx_plot, width = 10, height = 6, dpi = 200)
  output_filepath
}


#### added for triple rmsx, rmsd rmsf plot

plot_rmsd <- function(rmsd, Frame, RMSD) {
  rmsd_plot <- ggplot(rmsd, aes(x=Frame, y=RMSD)) +
    geom_line() +
    theme_minimal() +
    # coord_cartesian(xlim=c(0,max(rmsd$Frame))) # doesn't seem to do anything once it's there
    # xlim(, max(rmsd$Frame)) +
    labs(x = "")
  rmsd_plot
}




plot_rmsf <- function(rmsf_whole_traj, ResidueID, RMSF) {
  rmsf_plot <- ggplot(rmsf_whole_traj, aes(x=ResidueID, y=RMSF)) +
    geom_line() +
    theme_minimal() +
    coord_flip() +
    labs(x = "")
  rmsf_plot
}



plot_triple <- function(rmsx_plot, rmsd_plot, rmsf_plot) {
  plot_elements <- list(rmsd_plot, rmsx_plot, rmsf_plot, nullGrob())
  
  # Arrange the plots according to the layout matrix
  arranged_elements <- arrangeGrob(grobs = plot_elements, layout_matrix=
                                     rbind(c(4,4,1,1,1,1,1,1,4,4,4),
                                           c(4,4,1,1,1,1,1,1,4,4,4),
                                           c(4,2,2,2,2,2,2,2,3,3,3),
                                           c(4,2,2,2,2,2,2,2,3,3,3),
                                           c(4,2,2,2,2,2,2,2,3,3,3),
                                           c(4,2,2,2,2,2,2,2,3,3,3))) 
  
  # Return the arranged grob
  arranged_elements
  
}





# Main function to execute all steps
main <- function() {
  args <- parse_args()
  rmsx_raw <- read_and_summarize_csv(args$csv_path)
  for (id in unique(rmsx_raw$ChainID)) {
    rmsx_plot <- process_data_by_chain_id(rmsx_raw, id, args$csv_path, args$interpolate)
    save_plot(rmsx_plot, args$csv_path, id)
    
    # Ensure these data are properly loaded and variables are correctly defined
    rmsd_data <- read_csv(args$rmsd)  # Assuming args$rmsd is a path
    rmsf_data <- read_csv(args$rmsf)  # Assuming args$rmsf is a path
    rmsd_data
    rmsf_data
    rmsd_plot <- plot_rmsd(rmsd_data, rmsd_data$Frame, rmsd_data$RMSD)
    rmsf_plot <- plot_rmsf(rmsf_data, rmsf_data$ResidueID, rmsf_data$RMSF)
    triple_plot <- plot_triple(rmsx_plot, rmsd_plot, rmsf_plot)
    save_plot(triple_plot, args$csv_path, id)
    # Define the output file path
    output_file <- "./triple_plot.png"
    
    # Save the plot
    png(output_file, width = 800, height = 600)         # Open a PNG device
    grid::grid.draw(triple_plot)                        # Draw the plot
    dev.off()                                            # Close the device
    
    # Print the absolute path where the plot was saved
    cat("Plot saved to:", normalizePath(output_file), "\n")
    
    triple_plot
  }
}

main()




