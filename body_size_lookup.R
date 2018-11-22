#############################################################################
### Fill in body size information for all of the species in the community ###
#############################################################################

### Load body size data
all_bs_dat            <- read.csv("data/bs_data", sep = "")
if (load_matched_bs == TRUE & file.exists("saved_matching/needed_bs")) {
previously_matched_bs <- read.csv("saved_matching/needed_bs", sep = "")
} else {
previously_matched_bs <- NULL
}

### Data.frame for species body sizes
community_spec_bs <- data.frame(
  species          = c(
    needed_tips[["phylo_name"]]
  , have_spec_titer_surv[["phylo_name"]]
  , have_spec_bite[["phylo_name"]]
  , have_spec_detect[["phylo_name"]]
  , samp_data_com[which_rows_missing, ][["scientific_name"]])
, body_size        = 0)

## Determine whether the scientific name was matched or not
community_spec_bs <- transform(community_spec_bs
  , phylo_name_found = c(rep(1, nrow(community_spec_bs) -
      length(samp_data_com[which_rows_missing, ][["scientific_name"]]))
    , rep(0, length(samp_data_com[which_rows_missing, ][["scientific_name"]]))))

community_spec_bs <- community_spec_bs[!duplicated(community_spec_bs$species), ]

temp_spec_bs <- data.frame(
   phylo_name           = c(
     needed_tips[["phylo_name"]]           
   , have_spec_titer_surv[["phylo_name"]]  
   , have_spec_bite[["phylo_name"]]  
   , have_spec_detect[["phylo_name"]])
 , ebird_name           = c(
     as.character(needed_tips[["ebird_name"]]) 
   , as.character(have_spec_titer_surv[["ebird_name"]])
   , as.character(have_spec_bite[["ebird_name"]])
   , as.character(have_spec_detect[["ebird_name"]]))
 , body_size_ebird_name = 0
 , body_size_phylo_name = 0)

temp_spec_bs <- temp_spec_bs[!duplicated(temp_spec_bs[["ebird_name"]]), ]
temp_spec_bs <- temp_spec_bs[!duplicated(temp_spec_bs[["phylo_name"]]), ]

## First-level match, check both ebird names and phylo names, then fill in missing
temp_spec_bs[["body_size_ebird_name"]] <- all_bs_dat[
  match(as.character(temp_spec_bs$ebird_name)
    , as.character(all_bs_dat[["Scientific_name"]])), ][["body_size"]]

temp_spec_bs[["body_size_phylo_name"]] <- all_bs_dat[
  match(as.character(temp_spec_bs[["phylo_name"]])
    , as.character(all_bs_dat[["Scientific_name"]])), ][["body_size"]]

## fill in the missing values from one or the other lookup
temp_spec_bs[is.na(temp_spec_bs[["body_size_phylo_name"]]), 4] <- 
  temp_spec_bs[is.na(temp_spec_bs[["body_size_phylo_name"]]), 3]

temp_spec_bs[is.na(temp_spec_bs[["body_size_ebird_name"]]), 3] <- 
  temp_spec_bs[is.na(temp_spec_bs[["body_size_ebird_name"]]), 4]

## First step is to check previously saved matches to avoid lookup
missing_bs_dat <- which(is.na(temp_spec_bs[["body_size_ebird_name"]]))
unknown_names <- list(
  as.character(temp_spec_bs[is.na(temp_spec_bs[["body_size_ebird_name"]]), ][["phylo_name"]])
, as.character(temp_spec_bs[is.na(temp_spec_bs[["body_size_ebird_name"]]), ][["ebird_name"]]))

for (z in 1:length(unknown_names[[1]])) {
  
  prev_row <- which(
    as.character(previously_matched_bs[["species"]]) == unknown_names[[1]][z] |
    as.character(previously_matched_bs[["species"]]) == unknown_names[[2]][z])

  ## If the name cant be found at all...
  if (length(prev_row) == 0) prev_row <- NA
  
  ## if a matched name was found
  if (!is.na(prev_row)) {
  temp_spec_bs[missing_bs_dat[z], 3] <- previously_matched_bs[prev_row, ][["body_size"]]
  temp_spec_bs[missing_bs_dat[z], 4] <- previously_matched_bs[prev_row, ][["body_size"]]
  } 

}

### IUCN lookup for the missing species
missing_bs_dat <- which(is.na(temp_spec_bs[["body_size_ebird_name"]]))
unknown_names  <- as.character(temp_spec_bs[is.na(temp_spec_bs[["body_size_ebird_name"]]), ][["phylo_name"]])

## Each step in the next three loops overwrites unknown_names, pairing it down 
if (length(unknown_names > 1)) {

for (z in 1:length(unknown_names)) {

## here is where I need to search for a species's scientific name
 temp_nam <- unknown_names[z]
 
 ## IUCN synonym lookup
 if (!is.null(IUCN_REDLIST_KEY)) {
 temp_search <- rl_synonyms(name = temp_nam
  , key   = IUCN_REDLIST_KEY
  , parse = TRUE);  print(temp_search) ## print to keep track of what is going on
 
## print to keep track of what is going on
 print(temp_search) 
 
## or, just store an object with the same structure filled with NULL
 } else {
   temp_search <- list(
     name   = temp_nam
  ,  count  = 0
  ,  result = NULL)
 }

## set logicals/counters to initiate while loop
 temp_tip <- numeric(0)
 pp <- 1
 
## Make sure the IUCN lookup found an option
 if (class(nrow(temp_search[["result"]])) == "integer") {
 
## try each synonym in the phylogeny to find the correct one
 while (length(temp_tip) < 1) {
   
 new_nam <- temp_search[["result"]][pp, ][["synonym"]]
 new_nam <- strsplit(new_nam, split = " ")[[1]]
 new_nam <- paste(new_nam[1], new_nam[2], sep = "_")
 temp_tip <- which(all_bs_dat[["Scientific_name"]] == new_nam)
 
## counter to try next row of the synonym list
 pp <- pp + 1
 
 ## exit loop if search has made it to the end of the synonym search and there are no matches
 if (pp > nrow(temp_search$result)) break 
 
 }
   
 }
 
 ## however, sometimes it will fail. If so try a catalogue of life search
 if (length(temp_tip) < 1) {
   
 temp_search <- col_search(temp_nam, response = "terse")
 
 ## new counter
 pp <- 1
 
 if (nrow(temp_search[[1]]) >= 1) {
   
 ## extract the alternative names at the species level
 temp_search <- temp_search[[1]] %>% 
   dplyr::select(name, acc_name, rank, status) %>%
   filter(rank == "species")
 
## try each row in this lookup to find synonyms, continue until a synonym is found
 temp_tip <- numeric(0)
 while (length(temp_tip) < 1) {

 new_id   <- temp_search[pp, ] %>% dplyr::select(acc_name)
 new_nam  <- new_id[["acc_name"]]
 new_nam  <- strsplit(new_nam, split = " ")[[1]]
 new_nam  <- paste(new_nam[1], new_nam[2], sep = "_")
 temp_tip <- which(MyTree[["tip.label"]] == new_nam)
 
 ## counter to try next row of the synonym list
 pp <- pp + 1 
 
 ## exit loop if search has made it to the end of the synonym search and there are no matches
 if (pp > nrow(temp_search)) break 
 
 }
  
 ## end of else option for searching each row of id # pp
 } 

 ## end of Catalogue of Life lookup option  
 } 
 
  ## If the name in the phylogeny still isnt found
 if (length(temp_tip) < 1) {
   
 ## no names found for this choice in two databases, manual entry

 ## If running in an automated way from the command line, can set an option to ignore those species
 ## for which the automated procedure didnt find a fit
  if (skip_manual_input == FALSE) {
 
 ## let the user know they need to enter a value manually if the code is running in the background  
 beep(sound = 2, expr = NULL) 
  
 new_nam <- readline(prompt = "No synonyms found for this species, Manually Enter in the form Genus species (no _): ")
 new_nam <- strsplit(new_nam, split = " ")[[1]]
 new_nam <- paste(new_nam[1], new_nam[2], sep = "_")
 temp_tip <- which(MyTree[["tip.label"]] == new_nam)
 
 while (length(temp_tip) < 1) {
   
  new_nam <- readline(prompt = "Incorrect Name, Try Again. Type 'skip' to skip over this species: ")
  if (new_nam == "skip") break
  
  new_nam <- strsplit(new_nam, split = " ")[[1]]
  new_nam <- paste(new_nam[1], new_nam[2], sep = "_")
  temp_tip <- which(MyTree[["tip.label"]] == new_nam)
 }
   } 
    
 ## If a name was finally found, replace it
 if (length(temp_tip) != 0) {
   
 temp_spec_bs[missing_bs_dat[z], 3] <- all_bs_dat[temp_tip, ][["body_size"]]
 temp_spec_bs[missing_bs_dat[z], 4] <- all_bs_dat[temp_tip, ][["body_size"]] 
 
 ## Otherwise, stick in an NA for a name that wasn't found. These will be replaced later
 ## by Unknown, and estimated with body size (if it can be found) but no phylogenetic
 ## relationship
 } else {
 temp_spec_bs[missing_bs_dat[z], 3] <- NA
 temp_spec_bs[missing_bs_dat[z], 4] <- NA  
 }
 
 }
 
 ## IUCN will complain and then lock you out if you try too many searches in too short of a time.
 print("Waiting 10 seconds to escape being flagged as a bot"); Sys.sleep(10)
   
}
 
### Manually imput searching common names 
missing_bs_dat <- which(is.na(temp_spec_bs$body_size_ebird_name))
unknown_names  <- as.character(temp_spec_bs[is.na(temp_spec_bs$body_size_ebird_name), ]$phylo_name)

}

if (length(unknown_names > 1)) {

for (z in 1:length(unknown_names)) {
  
print(unknown_names[z])
  
new_nam  <- readline(prompt = "No synonyms found for this species, Manually Enter Common Name (spaces, no _): ")
temp_tip <- which(all_bs_dat$Common_name == new_nam)

 while (length(temp_tip) < 1) {
  new_nam  <- readline(prompt = "Incorrect Name, Try Again. Type 'skip' to skip over this speceies (spaces, no _): ")
  if (new_nam == "skip") break
  temp_tip <- which(all_bs_dat$Common_name == new_nam)
 }

 if (new_nam != "skip") {
 temp_spec_bs[missing_bs_dat[z], 3] <- all_bs_dat[temp_tip, ][["body_size"]]
 temp_spec_bs[missing_bs_dat[z], 4] <- all_bs_dat[temp_tip, ][["body_size"]] 
 } else {
 temp_spec_bs[missing_bs_dat[z], 3] <- NA
 temp_spec_bs[missing_bs_dat[z], 4] <- NA  
 }

}

### For all of those species with NAs still manually ask for body size
missing_bs_dat <- which(is.na(temp_spec_bs$body_size_ebird_name))
unknown_names  <- as.character(temp_spec_bs[is.na(temp_spec_bs$body_size_ebird_name), ]$phylo_name)

}

if (length(unknown_names > 1)) {

for (z in 1:length(unknown_names)) {
  
  print(unknown_names[z])
  
manually_entered_bs <- readline(prompt = "No species found. Enter body size (a reasonable option may be the mean of all species in this species' genus): ")


 temp_spec_bs[missing_bs_dat[z], 3] <- manually_entered_bs
 temp_spec_bs[missing_bs_dat[z], 4] <- manually_entered_bs
   
}
  
}

## Manually fill in body sizes for species without a matched name
rows_missing_names <- which(community_spec_bs[["phylo_name_found"]] == 0)
which_missing_name  <- community_spec_bs[rows_missing_names, ][["species"]]
if (length(which_missing_name) > 0) {
  
for (z in seq_along(which_missing_name)) {
  
    print(which_missing_name[z])
  
manually_entered_bs <- readline(prompt = "No species found. Enter body size (a reasonable option may be the mean of all species in this species' genus): ")

 community_spec_bs[rows_missing_names[z], 2] <- manually_entered_bs
 
}

}

### Obtain the consensus body sizes. Both columns should be the same, but take mean just in case...
community_spec_bs <- transform(community_spec_bs
  , body_size = as.numeric(c(
    rowMeans(
    data.frame(as.numeric(temp_spec_bs[["body_size_ebird_name"]])
  , as.numeric(temp_spec_bs[["body_size_phylo_name"]])))
  , community_spec_bs[community_spec_bs[["phylo_name_found"]] == 0, ][["body_size"]])))

## Export the matches for future data sets
if (overright_matched_bs == TRUE) {
  
  temp_full_list <- rbind(previously_matched_bs, community_spec_bs)
  temp_full_list <- temp_full_list[!duplicated(temp_full_list[, 1]), ]
  
  write.table(community_spec_bs, file = "saved_matching/needed_bs")
}

## convert to kg from g for lm 
community_spec_bs <- community_spec_bs %>%
  mutate(body_size_s = body_size / 1000)

### Add matched body size data into my data frames for running the models
temp_bs <- match(titercurves_reduced[["Scientific_Name"]], community_spec_bs[["species"]])
titercurves_reduced <- titercurves_reduced %>% 
  transform(body_size = community_spec_bs[temp_bs, ][["body_size"]]) %>%
  mutate(body_size_s  = body_size / 1000)

temp_bs <- match(survival_reduced[["Scientific_Name"]], community_spec_bs[["species"]])
survival_reduced <- survival_reduced %>% 
  transform(body_size = community_spec_bs[temp_bs, ][["body_size"]]) %>%
  mutate(body_size_s  = body_size / 1000)

temp_bs <- match(have_species_detect_full[["Scientific_Name"]], community_spec_bs[["species"]])
have_species_detect_full <- have_species_detect_full %>% 
  transform(body_size = community_spec_bs[temp_bs, ][["body_size"]]) %>%
  mutate(body_size_s  = body_size / 1000)

## Note---At this point species:
  ## (1) With a matched phylo name have a matched body size if it could be found, or a manually entered body size
  ## (2) Without a phylogeny entry but with a body size will be missing from "needed tips" but will have a scientific name
    ## recorded in community_spec_bs with a manually entered body size, which will be used in estimation
  ## (3) Without a phylogeny entry and no body size will be listed as Missing and estimated at the grand mean (this
    ## shouldn't occur using the above scripts but could happen if manual entry is ignored)
