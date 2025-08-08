# Minimal dependencies for plotting
pkgs <- c("ggplot2","viridis","dplyr","tidyr","stringr","readr","gridExtra")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos="https://cloud.r-project.org", quiet=TRUE)
cat("R plotting deps ready.\n")
