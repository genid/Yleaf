#!/usr/bin/env Rscript


packages <- c("mice", "plyr", "foreach")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}else {
  message("All required R packages are already installed, cool!")
}
