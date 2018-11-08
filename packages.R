#########################
### Required packages ###
#########################

needed_packages <- c(
  "lme4"
, "ape"
, "phytools"
, "ggplot2"
, "reshape2"
, "broom"
, "rstan"
, "tidyverse"
, "gridExtra"
, "taxize"
, "myTAI"
, "rredlist"
, "beepr"
, "data.table"
, "MASS"
, "PBSmapping"
, "maptools"
, "geosphere"
, "GISTools"
, "sf"
, "boot"
, "mgcv"
, "geosphere"
, "nlmrt"
, "shinystan"
, "rgdal"
, "tmap")

## Check if the packages are installed. If they are not install them, then load them
if (length(setdiff(needed_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(needed_packages, rownames(installed.packages())))  
}

lapply(needed_packages, require, character.only = TRUE)
