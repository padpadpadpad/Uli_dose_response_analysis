# install all necessary packages 

# package list in scripts/analysis.R
pkgs <- c('tidyr',
          'dplyr',
          'ggplot2',
          'nlme',
          'readxl',
          'viridis')

# check if they are already installed
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

# install packages
if(length(new_pkgs)) install.packages(new_pkgs)

# install rstan separately as needs to be compiled from source
install.packages('rstan', type = 'source') 