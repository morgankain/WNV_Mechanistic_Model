######################################
## Prepare temperature data for use ##
######################################

## Separate out the first column of the temperature data according to the details given in the 
 ## readme text on ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/

## 1-2: state code
## 3-5: county number
## 6-7: precip or temp code
## 8-11: year

temp_file_name <- paste("saved_output", paste(which_state, "county_temp_data", sep = "_"), sep = "/")

if (file.exists(paste(temp_file_name, ".Rds", sep = "")) & load_temp_data == TRUE) {
  
  county_temp_data <- readRDS(paste(temp_file_name, ".Rds", sep = ""))
  
} else {

## Actual temp data from: ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/
county_temp_data <- read.table("data/county_data/Texas_County_Temperature_Data", quote="\"", comment.char="")
## ID codes for each county from: https://www.dshs.texas.gov/chs/info/info_txco.shtm
county_ID        <- read.csv("data/county_data/Texas_FIPS_Codes.csv")
temp_split       <- strsplit(as.character(county_temp_data[, 1]), split = "")
temp_dat         <- data.frame(state = numeric(0), county = numeric(0), year = numeric(0))

## string split the first column into three useful columns
for (i in 1:length(temp_split)) {
temp_dat <- rbind(temp_dat
, data.frame(
   state  = as.numeric(paste(temp_split[[i]][1:2], collapse = ""))
 , county = as.numeric(paste(temp_split[[i]][3:5], collapse = ""))
 , year   = as.numeric(paste(temp_split[[i]][8:11], collapse = ""))
  )
)

if ((i/500) %% 1 == 0) { print(i / length(temp_split))}

}

## New data frame with temp data
county_temp_data <- cbind(temp_dat, county_temp_data[, -1])
names(county_temp_data)[4:15] <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

county_temp_data <- transform(county_temp_data, county = as.character(county))

## Convert Texas numeric counties to their names using their FIPS codes
county_ID_split <- strsplit(as.character(county_ID[, 2]), split = "")
temp_dat        <- data.frame(county_id = numeric(0))

for (i in 1:length(county_ID_split)) {
temp_dat <- rbind(temp_dat
, data.frame(
   county_id  = as.numeric(paste(county_ID_split[[i]][3:5], collapse = ""))
   )
  )
}

county_ID        <- cbind(county_ID[, 1], temp_dat)
names(county_ID) <- c("county_name", "county")
county_ID        <- transform(county_ID, county = as.character(county))

## combine the temp data with the county names
county_temp_data <- left_join(county_ID, county_temp_data, by = "county")

## remove unecessary columns and remove unused years
county_temp_data <- county_temp_data %>% 
  dplyr::select(-county, -state) %>% 
  filter(year > 1999 & year < 2019)

## change name to conform to the rest of the scripts
names(county_temp_data)[1]   <- "county"
## convert to long form
county_temp_data             <- melt(county_temp_data, c("county", "year"))
names(county_temp_data)[3:4] <- c("month", "temp_f")

## add temp in c 
county_temp_data <- transform(county_temp_data, temp = (temp_f - 32) * (5/9))

## round c to the nearest degree
county_temp_data <- transform(county_temp_data, temp = round(temp))
county_temp_data <- county_temp_data %>% dplyr::select(-temp_f)
county_temp_data <- transform(county_temp_data, month = as.character(month))

## Make sure month column doesn't have Xs
county_month_split <- strsplit(as.character(county_temp_data$month), split = "")
for (i in 1:length(county_month_split)) {
county_temp_data$month[i] <- as.numeric(paste(county_month_split[[i]][-1], collapse = ""))
}

## Save the output so it doesn't have to be run every time
saveRDS(county_temp_data, paste(temp_file_name, ".Rds", sep = ""))

} 
