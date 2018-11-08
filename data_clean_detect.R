#######
## species for which I have data, maximum detection distance
#######
have_species_detect_full <- read.csv("data/bird_detectability.csv", sep = ",")
have_species_detect_full <- have_species_detect_full %>% filter(!is.na(Detection_Distance))

names(have_species_detect_full)[4] <- "Scientific_Name"
