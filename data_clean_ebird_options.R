##################################################
### Options for loading and sorting ebird data ###
##################################################

restrict_location <- FALSE                 ## restrict by lat lon?
restrict_year     <- FALSE                 ## restrict year
  if (restrict_year == TRUE) {
year_range        <- c(2005, 2007)         ## subset the data to just look at dates from X-Y
                                           ## here 2005-2007 (years of Hamer et al. 2009 sampling)
  } else {
year_range        <- c("all", "all")
  }
  if (restrict_location == TRUE) {
lat_range         <- c(41.5667, 41.8333)   ## subset data for just a lat lon range 
                                           ## here 08' on each side of 41°42′N (in decimal form) (location of Hamer et al. 2009 sampling)
lon_range         <- c(87.6000, 87.6667)   ## here 08' on each side of 87°44′W (in decomal form) (location of Hamer et al. 2009 sampling)
  } else {
lat_range         <- c("all", "all")
lon_range         <- c("all", "all")
  }
file_name         <- "ebird_data_for_R/"   ## check the number of files in the ebird data folder
ebird_file_list   <- list.files(file_name) ## sequence over the number of ebrid_data_X files
  if (length(ebird_file_list) == 1) {      ## single or multiple independent chunks of ebird data
sing_file         <- TRUE
  } else {
sing_file         <- FALSE 
  }
