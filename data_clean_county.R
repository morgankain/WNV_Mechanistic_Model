#######################################################################
### Add county level covariates: habitat type, population size      ###
### Pipeline for a given county requires .shp data files            ### 
### Before running, the path needs to be set in top_level_script.R  ###
### to the folders with the shape files. This step can be skipped   ###
### to estimate R0 for a given community without fitting a spatio   ###
### temporal model                                                  ###
#######################################################################

county_file_name <- paste("saved_output", paste(which_state, "counties_summary", sep = "_"), sep = "/")

## If the county information has been collated previously and the user does not want to rerun code, load data
if (file.exists(county_file_name) & load_prev_county_dat == TRUE & fit_spatio_temporal == TRUE) {
  
  county_data <- read.table(county_file_name, sep = "")

} else if (fit_spatio_temporal == TRUE) {

sids <- readShapePoly(county_boundaries_path)

########
## Calculate the geographic centroid of each county
########

writeSpatialShape(sids, "sids")
cents     <- coordinates(sids)
cents     <- SpatialPointsDataFrame(coords=  cents, data = sids@data, proj4string = CRS("+proj=longlat +ellps=clrk66"))
writeSpatialShape(cents, "cents")
centroids        <- getSpPPolygonsLabptSlots(sids)
centroids        <- as.data.frame(centroids)
names(centroids) <- c("X", "Y")
centroids$county <- as.character(unlist(sids[[1]]))
Current_Counties <- centroids

########
## Vegetation/habitat type map layer, determined as the proportion of each county covered by type X
########

hab_layer_info <- importShapefile(county_region_type_path)
hab_layer      <- readShapePoly(county_region_type_path)

hab_overlap      <- raster::intersect(sids, hab_layer)
hab_overlap$area <- raster::area(hab_overlap) / 1000000
hab_by_county    <- aggregate(area ~ CNTY_NM + REGIONS, data = hab_overlap, FUN = sum)

hab_by_county_sum <- hab_by_county %>%
  group_by(CNTY_NM) %>%
  mutate(total_area = sum(area))

hab_by_county_sum <- transform(hab_by_county_sum,
 prop_hab_type = area / total_area)

all_hab_type <- hab_by_county_sum[, -c(3, 4)] %>%
  spread(REGIONS, prop_hab_type, fill = 0)

## stick results together
which_rows       <- match(Current_Counties[["county"]], all_hab_type[["CNTY_NM"]])
Current_Counties <- cbind(Current_Counties, all_hab_type[which_rows, -1])

########
## Population census of each county
########

county_data <- read.csv(county_pop_density_path)
county_data <- transform(county_data, County = as.character(County))

which_rows  <- match(Current_Counties$county, county_data$County)

county_data <- cbind(county_data[which_rows, ]
  , dplyr::select(Current_Counties, -county))

## Double check to make sure data is in the correct format
if (class(county_data$Density) != "numeric") {
  county_data$Density <- as.numeric(as.character(county_data$Density))
}
if (class(county_data$Density) != "character") {
  county_data$County  <- as.character(county_data$County)
}

rownames(county_data) <- NULL

## Small Texas county naming problem to match the ebird data (de Witt doesn't have a space)
county_data <- transform(county_data, County = as.character(County))
county_data[["County"]][which(county_data[["County"]] == "De Witt")] <- "DeWitt"

## Save file so that all of the above doesn't need to be run each time
write.table(county_data, file = county_file_name)

}
