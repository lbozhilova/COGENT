################################################
##### COGENT tutorial package installation #####
################################################

# List all necessary CRAN packages
requiredPackages <- c("devtools",
                      "ggplot2",
                      "ggthemes",
                      "igraph")

# Check whether they are already installed
notInstalled <- setdiff(requiredPackages, installed.packages()[,"Package"])

# Install missing packages
if (length(notInstalled) > 0)
  install.packages(notInstalled)

# Install COGENT if you don't already have it
if(!("COGENT" %in% installed.packages()[,"Package"])){
  require("devtools")
  install_github("lbozhilova/COGENT")
}

# Clean up
rm(notInstalled, requiredPackages); gc()
