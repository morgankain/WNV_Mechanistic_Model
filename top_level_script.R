###############################################################################################
### Code for: Kain and Bolker 2019: Predicting West Nile virus transmission in North        ###
### American bird communities using phylogenetic mixed effects models and eBird citizen     ###
### science data. For questions email: kainm@mcmaster.ca OR morganpkain@gmail.com           ###
###############################################################################################

###############################################################################################
### If you have retrieved this code from github: all data needed to run these models can be ###
### obtained in the supplementary material of Kain and Bolker 2019 at:                      ###
### URL HERE                                                                                ###
###############################################################################################

###############################################################################################
### This code covers: Steps to estimating WNV community transmission using phylogenetic     ###
### mixed effects models, eBird data and data from the literature on bird  and mosquito     ###
### physiological responses to WNV infection, mosquito biting preferences, bird             ###
### detectability and other data                                                            ###
###############################################################################################

## What state?
which_state <- "Texas"

##############
### Some controls
##############
skip_manual_input       <- FALSE            ## For automated shell executing of code avoid manual entry options (all scientific bird names won't get matched)
test_species_importance <- FALSE            ## Run model for species importance? (takes considerable time, so may want to turn off)

##############
### Packages and functions
##############
source("packages.R")                        ## required packages
source("ggplot_theme.R")                    ## theme for ggplots
source("phylo_setup.R")                     ## phylo_lmm, phylo_glmm and helper functions
source("phylo_prediction_tools.R")          ## functions for fitting models and estimating species' responses

##############
### Set up what uncertainty is propagated
##############
source("uncertainty_setup.R")

## Check if there is no uncertainty (no uncertainty reduces how estimates are passed)
no_uncer <- ifelse(sum((uncertainty_list == FALSE)) == 0, FALSE, TRUE)

##############
### Read in and clean previously gathered data for models
##############
source("phylo_clean.R")                     ## Phylogeny, set up to load previously saved consensus phylogeny
                                            ## User calculated phylogeny can be substituted with naming convention "consensus"
source("data_clean_titer.R")                ## Bird titer profiles
source("data_clean_survival.R")             ## Bird survival
source("data_clean_biting.R")               ## Mosquito biting preferences
source("data_clean_detect.R")               ## Bird detection scaling
load_temp_data          <- TRUE
source("data_clean_temperature.R")          ## State level county data, already subset to Texas. Full dataset available at: ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/
load_stan_bird_mos_res  <- TRUE             ## Load previously run model if available (TRUE) or run from scratch (FALSE)?
source("data_clean_bird_mos.R")             ## Bird to Mosquito transmission
load_stan_mos_bird_res  <- TRUE             ## Load previously run model if available (TRUE) or run from scratch (FALSE)?
source("data_clean_mos_bird.R")             ## Mosquito to Bird transmission
load_prev_county_dat    <- TRUE             ## Load previously run model if available (TRUE) or run from scratch (FALSE)?
fit_spatio_temporal     <- TRUE             ## Is data available to fit a spatio-temporal model?
                                            ## If TRUE, and paths supplied are wrong, an error will occur
                                            ## Path to county boundaries
county_boundaries_path  <- "data/county_data/txdot_counties/txdot_2015_county_detailed_tx.shp"
                                            ## Path to county regions
county_region_type_path <- "data/county_data/naturalregions_tx/Natural_Regions_Maj.shp"
                                            ## Path to county pop density
county_pop_density_path <- "data/county_data/Texas_pop_density.csv"
source("data_clean_county.R")               ## Gather and organize county data

##############
### Calculate biting preferences for the species from Hamer et al. 2009 Using a Bayesian approach.
### Predict biting preferences for both current and new species
### Needed prior to setting up and running phylo methods for biting preference
##############
load_stan_bite_pref_res <- TRUE             ## load previously run model or run?
source("biting_pref_stan_model.R")          ## Stan biting model

##############
### Read in and clean ebird and phylogeny.
### Match scientific names from ebird data and body size data to phylogeny
##############
source("data_clean_ebird_options.R")        ## options for loading and sorting ebird data
load_ebird_dat                <- TRUE       ## Load previously run model if available (TRUE) or run from scratch (FALSE)?
source("data_clean_ebird.R")                ## clean and gather csv files downloaded from ebird
load_matched_phylo_names      <- TRUE       ## Load previously matched scientific names?
overright_matched_phylo_names <- FALSE      ## Overwrite saved csv with new matched scientific names?
                                            ## IUCN redlist key. If no key is available change this to NA
IUCN_REDLIST_KEY <- "0a2af13fab231679376c105ae0591f5050613cf6b7ea9faeea2224d478e0218f"
source("match_avian_scinames_IUCN.R")       ## Make sure the species names in this ebird data frame match the phylogeny data frame

###############
### Fill in body size info
###############
load_matched_bs               <- TRUE       ## Load previously searched body sizes?
overright_matched_bs          <- FALSE      ## Overwrite saved csv with matched phylo names?
                                            ## Turn on whenever a new data set is loaded
source("body_size_lookup.R")

##############
### Predict responses for species outside of the data set
##############
used_saved_responses_titer    <- TRUE       ## load previous responses for titer if possible?
                                            ## (given a previous saved result for the given uncertainty)
used_saved_responses_survival <- TRUE       ## load previous responses for survival if possible?
used_saved_responses_bite     <- TRUE       ## load previous responses for biting preferences if possible?
used_saved_responses_detect   <- TRUE       ## load previous responses for species detectability if possible?
source("new_responses_titer_surv.R")
source("new_responses_prop_bite_pref.R")
source("new_responses_detectability.R")

## Clear up some RAM
gc()

##############
### Calculate community competence, at both the species (unweighted by proportion in community) and community level (weighted by proportion in the community)
##############
source("species_competence.R")              ## Physiological competence of each species
source("community_competence.R")            ## Using physiological competence, biting preference etc. to calculate community R0

##############
### Free up some ram
##############

if (exists("comm_comp_summary_m")) {
  rm(comm_comp_summary_m)
}

## Clear up some RAM
gc()

##############
### Check correlations/dilution effect assumptions and save results
##############

source("summary_stats_cleaned.R")           ## Statistics on the results
source("saving_script.R")                   ## Save results to disk
rmarkdown::render('summary_output.Rmd')     ## Produce a pdf summary of results

###########################################################
### Exploration of results outside of the main analyses ###
###########################################################
# source("summary_stats_compare_cleaned.R")            ## determine magnitude of effect of uncertainty (after all runs complete)
# source("how_many_complete_lists.R")                  ## examine effects of resampling on species coverage and error in species proportions
