##########################################################################################
### Generic script to take ebird data, downloaded in the form of csv files, and subset ###
### them into appropriate dates/species etc.                                           ###
### Requires a folder in the main directory called ebird_data with each file in this   ###
### folder named ebird_data_X, where X is a sequence by 1 starting at 1.               ###
### This will be set up automatically if the ebird_data_clean.sh file is used          ###
##########################################################################################

################################################
### Script to be used for multiple locations ###
################################################

## check for a previously saved summarized file.
## check for appropriate month, year, county level summarization
temp_opts       <- paste(restrict_location, restrict_year, sep = "_")
ebird_file_name <- paste("saved_output", paste(which_state, "summarized_ebird_data", temp_opts, sep = "_"), sep = "/")
ebird_meta_name <- paste("saved_output", paste(which_state, "ebird_metadata", temp_opts, sep = "_"), sep = "/")

if (file.exists(ebird_file_name) & file.exists(ebird_meta_name) & load_ebird_dat == TRUE) {
  
  samp_data_com       <- read.csv(ebird_file_name, sep = "")
  effort_metadata_com <- read.csv(ebird_meta_name, sep = "")  
  
} else {

for (i in seq_along(ebird_file_list)) {

## Print progress
print(i/length(ebird_file_list))

## Read in data
samp_data <- read.delim(
    paste(file_name,  ebird_file_list[i], sep = "")
  , quote = ""
  , head = FALSE
  , row.names = NULL
  , stringsAsFactors = FALSE)
  
## If the first row is a header row, set it as the names of the columns and save the column names.
## Otherwise give the column names to the file from the stored column names
if (samp_data[1, 1] == "COMMON NAME") {
  ebird_col_names     <- samp_data[1, ]
  samp_data           <- samp_data[-1, ]
  colnames(samp_data) <- ebird_col_names
} else {
  colnames(samp_data) <- ebird_col_names
}
  
## break down the date column into their part and make three new columns
dates_vec <- strsplit(as.character(as.Date(samp_data[['OBSERVATION DATE']])), split = "-")
years_vec <- numeric(length(dates_vec))
month_vec <- numeric(length(dates_vec))
for (j in 1:length(years_vec)) {
  years_vec[j] <- dates_vec[[j]][1]
  month_vec[j] <- dates_vec[[j]][2]
}

## After a transform statement, spaces in column names have been replaced with dots "."
samp_data <- transform(samp_data
  , YEAR  = as.numeric(years_vec)
  , MONTH = as.numeric(month_vec)
  , DAY   = as.numeric(format(as.Date(samp_data[['OBSERVATION DATE']]), "%j")))

rm(dates_vec); rm(years_vec); rm(month_vec)

## subset the data or simply rename it
if (restrict_year == TRUE) {
  samp_data <- samp_data %>% filter(YEAR >= year_range[1], samp_data$YEAR <= year_range[2])
} 

## clean up the data frame for use (remove those observations that do not have a count)
samp_data <- samp_data %>% 
  transform(OBSERVATION.COUNT = as.numeric(OBSERVATION.COUNT)) %>%
  filter(!is.na(OBSERVATION.COUNT))

## select the area
if (restrict_location == TRUE) {

 samp_data <- samp_data %>% 
   transform(LONGITUDE = LONGITUDE*-1) %>%
   filter(
     LATITUDE >=  lat_range[1]
   , LATITUDE <=  lat_range[2]
   , LONGITUDE >= lon_range[1]
   , LONGITUDE <= lon_range[2])

}

## Need properties of a data frame not tibble
samp_data <- as.data.frame(samp_data)

## remove all counts that are not for a single species (unknowns, hybrids etc.) (Minimal # of counts)
inc_rem <- grep("sp[.]", samp_data[['SCIENTIFIC.NAME']])
if (length(inc_rem) >= 1) {
  samp_data <- samp_data[-inc_rem, ]
}
inc_rem <- grep("hybrid", samp_data[['SCIENTIFIC.NAME']])
if (length(inc_rem) >= 1) {
  samp_data <- samp_data[-inc_rem, ]
}
inc_rem <- grep("/", samp_data[['SCIENTIFIC.NAME']])
if (length(inc_rem) >= 1) {
  samp_data <- samp_data[-inc_rem, ]
}
inc_rem <- grep(" x ", samp_data[['SCIENTIFIC.NAME']])
if (length(inc_rem) >= 1) {
  samp_data <- samp_data[-inc_rem, ]
}
inc_rem <- grep("type", samp_data[['SCIENTIFIC.NAME']])
if (length(inc_rem) >= 1) {
  samp_data <- samp_data[-inc_rem, ]
}

## add some identification to the data frame and save it 
samp_data <- transform(samp_data
  , lat_min  = lat_range[1]
  , lat_max  = lat_range[2]
  , lon_min  = lon_range[1]
  , lon_max  = lon_range[2]
  , date_min = year_range[1]
  , date_max = year_range[2])

## remove all non-complete lists (putting this at the end of cleaning instead of the beginning in 
 ## case a user would want to use non-complete lists as well, given that a number was given
  ## for each species count (not an X))
samp_data <- samp_data %>% filter(ALL.SPECIES.REPORTED == 1)

## combine each loaded csv
if (i == 1) {
 samp_data_com       <- samp_data
} else {
 samp_data_com       <- rbind(samp_data_com, samp_data)
}

}
  
###############
### Calculate Effort Meta-Data
###############
  
## Number of unique lists submitted and mean effort (distance and time)
  ## Combinations of County, Month, and Year that do not have these data recorded will be listed as NA
temp_effort_metadata <- samp_data_com %>% 
  group_by(COUNTY, MONTH, YEAR, SAMPLING.EVENT.IDENTIFIER) %>%
  summarize(mean_effort_dist = mean(as.numeric(EFFORT.DISTANCE.KM))
          , mean_effort_time = mean(as.numeric(DURATION.MINUTES)))

## Number of unique lists and the number of lists that record distance and time
  ## Count the proportions of lists in County, Month, and Year that are not NA (these are the lists that recorded
    ## distance and/or time data)
temp_effort_metadata <- temp_effort_metadata %>%
  group_by(COUNTY, MONTH, YEAR) %>%
  summarize(num_lists                = length(unique(SAMPLING.EVENT.IDENTIFIER))
          , prop_report_effort_dist  = length(which(!is.na(mean_effort_dist))) / length(mean_effort_dist)
          , prop_report_effort_time  = length(which(!is.na(mean_effort_time))) / length(mean_effort_time))

## Subset out the checklists that have distance data, and check the average distance effort
effort_metadata1 <- samp_data_com  %>% 
  filter(EFFORT.DISTANCE.KM != "") %>%
  group_by(COUNTY, MONTH, YEAR)    %>%
  summarize(mean_distance = mean(as.numeric(EFFORT.DISTANCE.KM), na.rm = TRUE))

## Subset out the checklists that have time data, and check the average time effort
effort_metadata2 <- samp_data_com  %>% 
  filter(DURATION.MINUTES != "")   %>%
  group_by(COUNTY, MONTH, YEAR)    %>%
  summarize(mean_time      = mean(as.numeric(DURATION.MINUTES), na.rm = TRUE))

## Merge meta-data
effort_metadata_m1 <- merge(effort_metadata1, temp_effort_metadata, by = c("COUNTY", "MONTH", "YEAR"), all = TRUE)
effort_metadata_m2 <- merge(effort_metadata2, effort_metadata_m1, by = c("COUNTY", "MONTH", "YEAR"), all = TRUE)

## Estimate total time assuming lists with unrecorded distance and time data are at the mean of those reported
effort_metadata <- transform(effort_metadata_m2
  , num_lists_report_dist = num_lists     * prop_report_effort_dist
  , num_lists_report_time = num_lists     * prop_report_effort_time
  , total_distance        = mean_distance * (1/prop_report_effort_dist)
  , total_time            = mean_time     * (1/prop_report_effort_time))

## rename to match bird list data frame 
effort_metadata_com <- effort_metadata; rm(effort_metadata)

## first save the complete samp_data_com for use in a separate simulation later
 ## (not needed for the main analysis)
samp_data_com_full <- samp_data_com
write.table(samp_data_com_full, file = paste(ebird_file_name, "FULL", sep = "_"))
rm(samp_data_com_full)

## After the loop over csv files is complete, aggregate the counts a second time
  ## needed because of the hard break in the csv files imposed, which can cut a
    ## single location/time into two csv files
samp_data_com <- samp_data_com %>% 
  group_by(COUNTY, MONTH, YEAR, SCIENTIFIC.NAME) %>%
  summarize(total_obs = sum(OBSERVATION.COUNT))

## Write the data to disk for use later
  write.table(samp_data_com, file = ebird_file_name)
  write.table(effort_metadata_com, file = ebird_meta_name)

}
  