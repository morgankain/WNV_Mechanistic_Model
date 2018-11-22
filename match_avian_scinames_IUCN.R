###################################################################################
### Automated script to match scientific names of the avian phylogeny and ebird ###
### using API token from ICUN. In the absence of an IUCN key users will have to ### 
### rely more heavily on manual input.                                          ### 
###################################################################################

## rename column headings from the ebird file
names(samp_data_com) <- c("county", "month", "year", "scientific_name", "ebird_prior")   

## adjust scientific names step 1 (change space to an underscore)
temp_nams     <- with(samp_data_com, tstrsplit(scientific_name, split = " "))
samp_data_com <- transform(samp_data_com, scientific_name = paste(temp_nams[[1]], temp_nams[[2]], sep = "_"))
rm(temp_nams)

#######
## tip numbers of the full avian phylogeny that correspond to the species that I need
#######

ebird_species_list <- unique(samp_data_com[['scientific_name']])

needed_tips <- data.frame(
  ebird_name = ebird_species_list
, phylo_name = numeric(length(ebird_species_list))
, tip_number = numeric(length(ebird_species_list)))

## Vector for the list of species that have different names
unknown_names <- numeric(0)

## First get all of the known names
pp <- 1
for (z in 1:nrow(needed_tips)) {
  
## Find which scientific name in the phylogeny matches with the name in the ebird species list
## Some of these will correspond to the same bird, but with different names, which will be retunred
## as NAs. The objective of this script is to resolve these naming differences
  temp_tip <- which(MyTree[["tip.label"]] == as.character(needed_tips[z, 1]))
  
  if (length(temp_tip) < 1) { 
    unknown_names[pp] <- z
    pp <- pp + 1
    
 needed_tips[z, 2] <- "na"
 needed_tips[z, 3] <- 0
  } else {
    
 needed_tips[z, 2] <- MyTree[["tip.label"]][temp_tip]
 needed_tips[z, 3] <- temp_tip
    
 }
  
}

## can get all of the unknown names by loading a previously saved data set
if (load_matched_phylo_names == TRUE & file.exists("saved_matching/matched_scinames")) {
  previously_matched_names <- read.csv("saved_matching/matched_scinames", sep = "")  
} else {
  previously_matched_names <- NULL
}

## second vector of unknown names
unknown_names2 <- numeric(0)
pp <- 1

## update the names with previously searched *and found* names
for (z in unknown_names) {
  
  ## find and try previously matched scientific name, one at a time
  prev_row <- which(as.character(previously_matched_names[["ebird_name"]]) == as.character(ebird_species_list[z]))
  
  ## If the name cant be found at all...
  if (length(prev_row) == 0) prev_row <- NA
  
  ## if a matched name was found
  if (!is.na(prev_row)) {
  found_tip         <- as.character(previously_matched_names[prev_row, ][["phylo_name"]])
  needed_tips[z, 2] <- found_tip 
  temp_tip          <- which(MyTree[["tip.label"]] == found_tip)
  needed_tips[z, 3] <- temp_tip
  
  ## if it wasn't, retain the record in a new vector
  } else {
    
  ## poor practice of extending a vector...
  unknown_names2[pp] <- z
  pp <- pp + 1
    
  }

}

## then check the remaining names
for (z in unknown_names2) {

## remove the "_" for automated search
 temp_nam <- strsplit(as.character(needed_tips[z, 1]), "_")[[1]]
 temp_nam <- paste(temp_nam[1], temp_nam[2], sep = " ")

## IUCN synonym lookup
 if (!is.null(IUCN_REDLIST_KEY)) {
 temp_search <- rl_synonyms(name = temp_nam
  , key   = IUCN_REDLIST_KEY
  , parse = TRUE)
 
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
 pp       <- 1
 
## Check if the IUCN lookup found an option
 if (class(nrow(temp_search[["result"]])) == "integer") {
 
## match each synonym returned by the IUCN with the phylogeny to find a match
 while (length(temp_tip) < 1) {
   
 ## return back to having an "_"
 new_nam  <- temp_search[["result"]][pp, ][["synonym"]]
 new_nam  <- strsplit(new_nam, split = " ")[[1]]
 new_nam  <- paste(new_nam[1], new_nam[2], sep = "_")
 temp_tip <- which(MyTree[["tip.label"]] == new_nam)
 
 ## counter to try next row of the synonym list
 pp <- pp + 1 
 
 ## exit loop if search has made it to the end of the synonym search and there are no matches
 if (pp > nrow(temp_search[["result"]])) break 
 
 }
   
 }

 ## If the IUCN search doesnt return a hit (or no key is available) try a Catalogue of Life search
 if (length(temp_tip) < 1) {
   
 ## Catalog of Life search
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
 needed_tips[z, 2] <- MyTree[["tip.label"]][temp_tip]
 needed_tips[z, 3] <- temp_tip
 ## Otherwise, stick in an NA for a name that wasn't found. These will be replaced later
 ## by Unknown, and estimated with body size (if it can be found) but no phylogenetic
 ## relationship
 } else {
 needed_tips[z, 2] <- NA
 needed_tips[z, 3] <- NA  
 }
 
 }
 
 ## IUCN will complain and then lock you out if you try too many searches in too short of a time.
 print("Waiting 10 seconds to escape being flagged as a bot"); Sys.sleep(10)
   
}

## save the unknown names for use again later, only adding the new names, if any have been added, removing duplicates
if (exists("previously_matched_names")) {
all_tips <- rbind(needed_tips, previously_matched_names)
} else {
all_tips <- needed_tips
}
all_tips <- all_tips[!duplicated(all_tips), ]

if (overright_matched_phylo_names == TRUE) {
  write.table(all_tips, file = "saved_matching/matched_scinames")
}

## Store the name from the phylogeny alongside the name from the ebird file
names(all_tips)[1] <- "scientific_name" 
samp_data_com      <- left_join(samp_data_com, all_tips, by = "scientific_name")

## Convert all NA values in the phylo_name column to be called "Missing", which are species 
  ## that will not obtain phylogenetically informed predictions. Later, this species will
    ## still be searched for by their scientific_name for body sizes
which_rows_missing <- which(is.na(samp_data_com[["phylo_name"]]))

if (length(which_rows_missing) != 0) {
  samp_data_com[which_rows_missing, ][["phylo_name"]] <- "Missing"
}

## Store which species were not found (for body size lookup), then remove the NA birds from 
  ## needed_tips because needed_tips will be used later to obtain phylogenetically informed estimates
which_rows_no_phylo    <- which(is.na(needed_tips[["phylo_name"]]))

## Record these species names, because they will still be searched later for body size
which_species_no_phylo <- needed_tips[which_rows_no_phylo, ]

if (length(which_rows_missing) != 0) {
needed_tips <- needed_tips[-which_rows_no_phylo, ]
}

#######
## Clean out clearly absurd observations.
## Not dynamic, use discretion
## User chosen species should be added here for other locations
#######

if (which_state == "Texas") {
  
  samp_data_com <- samp_data_com %>%
    filter(
       phylo_name != "Struthio_camelus" |
       phylo_name != "Dromaius_novaehollandiae" |
       is.na(phylo_name)
      )
  
}

#######
## species for which I have data, titer and survival
#######

## Ebird and phylo names are the same for these species
have_spec_titer_surv <- data.frame(
  ebird_name = needed_species[["Scientific.Name"]]
, phylo_name = needed_species[["Scientific.Name"]]
, tip_number = numeric(nrow(needed_species)))

for (z in 1:nrow(have_spec_titer_surv)) {
  
 temp_tip                   <- which(MyTree[["tip.label"]] == as.character(have_spec_titer_surv[["phylo_name"]])[z])
 have_spec_titer_surv[z, 3] <- temp_tip
    
  }

## find which species in the data set I do not have data for
missing_spec_titer_surv <- which(is.na(match(needed_tips[["phylo_name"]], have_spec_titer_surv[["phylo_name"]])))

have_spec_titer_surv <- transform(have_spec_titer_surv
  , ebird_name = as.character(ebird_name)
  , phylo_name = as.character(phylo_name))

#######
## species for which I have data, biting preference
#######

## Ebird and phylo names are the same for these species
have_spec_bite <- data.frame(
  ebird_name = as.character(bite_pref_res[["Scientific_Name"]])
, phylo_name = as.character(bite_pref_res[["Scientific_Name"]])
, tip_number = numeric(nrow(bite_pref_res)))

for (z in 1:nrow(have_spec_bite)) {
  
 temp_tip             <- which(MyTree[["tip.label"]] == as.character(have_spec_bite[["phylo_name"]])[z])
 have_spec_bite[z, 3] <- temp_tip
    
}

## find which species in the data set I do not have data for
missing_spec_bite <- which(is.na(match(needed_tips[["phylo_name"]], have_spec_bite[["phylo_name"]])))

have_spec_bite <- transform(have_spec_bite
  , ebird_name = as.character(ebird_name)
  , phylo_name = as.character(phylo_name))

#######
## species for which I have data, detection
#######

## phylo names were matched and exported as part of gathering data
have_spec_detect <- data.frame(
  ebird_name = as.character(have_species_detect_full[["SCIENTIFIC.NAME"]])
, phylo_name = as.character(have_species_detect_full[["Scientific_Name"]])
, tip_number = numeric(nrow(have_species_detect_full)))

for (z in 1:nrow(have_spec_detect)) {
  
 temp_tip               <- which(MyTree[["tip.label"]] == as.character(have_spec_detect[["phylo_name"]])[z])
 have_spec_detect[z, 3] <- temp_tip
    
}

## find which species in the data set I do not have data for
missing_spec_detect <- which(is.na(match(needed_tips[["phylo_name"]], have_spec_detect[["phylo_name"]])))

have_spec_detect <- transform(have_spec_detect
  , ebird_name = as.character(ebird_name)
  , phylo_name = as.character(phylo_name))

## just the vector of species names
spec_est_vec        <- unique(samp_data_com[["scientific_name"]])

  ## check if multiple communities or a single community will be analyzed
  sing_community <- FALSE

  ## Aggregate to a single list for efficient passing to a function later
specs_to_est_titer_surv <- list(
  needed_tips  = needed_tips[["tip_number"]]
, have_tips    = have_spec_titer_surv[["tip_number"]]
, missing_spec = missing_spec_titer_surv)

specs_to_est_bite <- list(
  needed_tips  = needed_tips[["tip_number"]]
, have_tips    = have_spec_bite[["tip_number"]]
, missing_spec = missing_spec_bite)

specs_to_est_detect <- list(
  needed_tips  = needed_tips[["tip_number"]]
, have_tips    = have_spec_detect[["tip_number"]]
, missing_spec = missing_spec_detect)
