######################################################################
### Calculate some summary stats, cleaned for automated script run ###
######################################################################

#######
## A parameter to adjust to designate if all spatial and temporal aggregations of the bird
## community and temperature. If only one of these aggregations is desired, scroll to the
## bottom of the script and run manually. Setting this parameter to TRUE increases the
## computation time by quite a lot
#######
run_spat_temp_aggregated <- FALSE

#######
## First, test some assumptions of the dilution effect
## Second, fit a spatio-temporal model
## Third, calculate the importance of species in a few different ways
#######

#######
## (1) Dilution effect assumption
#######

host_comp_summary_m <- melt(host_comp_summary, c("day", "log_dose", "species", "model", "outcome"))

host_comp_summary_p <- host_comp_summary_m %>% 
  filter(outcome == "realized_transmission" | outcome == "biting_pref") %>%
  spread(outcome, value); rm(host_comp_summary_m)

## Sort for testing dilution effect assumptions
species_to_est_data_p <- samp_data_com %>% 
  group_by(county, month, year) %>%
  summarize(tot_obs = sum(ebird_prior))

species_to_est_data_p <- merge(samp_data_com, species_to_est_data_p)

## Scale counts given detectability
if (no_uncer == TRUE) {

comm_detect_medians <- cbind(
   comm_detect_est[, c(1:3)]
  , data.frame(detectability = comm_detect_est[, -c(1:3)]))

} else {
  
comm_detect_medians <- cbind(
   comm_detect_est[, c(1:3)]
  , data.frame(detectability = apply(comm_detect_est[, -c(1:3)], 1, median)))
  
}

names(comm_detect_medians)[1] <- "phylo_name"
 
species_to_est_data_p <- left_join(
  species_to_est_data_p
  , comm_detect_medians[, c(1, 4)])

species_to_est_data_p <- species_to_est_data_p %>%
  group_by(county, month, year) %>%
  mutate(ebird_prior_adj = (1 / (detectability / max(detectability))) * ebird_prior)

species_to_est_data_p <- species_to_est_data_p %>%
  group_by(county, month, year) %>%
  mutate(tot_obs_adj = sum(ebird_prior_adj))

## Relative proportion of each bird in each community
species_to_est_data_p <- species_to_est_data_p %>%
  mutate(prop_obs_adj = ebird_prior_adj / tot_obs_adj)

## Merge in effort to be able to subset by # of lists
species_to_est_data_p <- left_join(species_to_est_data_p
  , effort_metadata_com[, c(1, 2, 3, 6)])

## reduce to communities with > 100 lists
species_to_est_data_p2 <- species_to_est_data_p %>% 
  filter(num_lists >= 80)

## median abundance of each species for complete and reduced ebird data sets
species_to_est_data_p <- species_to_est_data_p %>%
  group_by(scientific_name) %>%
  summarize(med_obs = median(prop_obs_adj))

species_to_est_data_p2 <- species_to_est_data_p2 %>%
  group_by(scientific_name) %>%
  summarize(med_obs = median(prop_obs_adj))

names(species_to_est_data_p)[1]  <- "species"
names(species_to_est_data_p2)[1] <- "species"

## combine adjusted communities with estimated species physiological competence
 ## host_comp_summary_p1 from check_species_importance.R
dilut_assump_full <- left_join(host_comp_summary_p2[, c(1, 2, 4)]
  , species_to_est_data_p)

dilut_assump_red <- left_join(host_comp_summary_p2[, c(1, 2, 4)]
  , species_to_est_data_p2)

## Lastly, add in temperature, which must be controlled for as it will covary strongly with
 ## species richness
comm_comp_summary_p2 <- transform(comm_comp_summary_p2, month = as.character(month))
comm_comp_summary_p2_well_observed <- transform(comm_comp_summary_p2_well_observed, month = as.character(month))

if (no_uncer == TRUE) {

dilut_eff_lm_red           <- lm(med_comp ~ log(med_obs)
  , data = dilut_assump_red)

dilut_eff_pred_lm_red      <- lm(med_comp ~ log(meannum_spec) + temp
  , data = comm_comp_summary_p2_well_observed)

} else {
  
dilut_eff_lm_red      <- lm(med_comp ~ log(med_obs)
  , weights = 1 / var_comp
  , data = dilut_assump_red)

dilut_eff_pred_lm_red <- lm(med_comp ~ log(meannum_spec) + temp 
  , weights = 1 / var_comp
  , data = comm_comp_summary_p2_well_observed)

dilut_eff_pred_lm <- lm(med_comp ~ log(meannum_spec) + temp 
  , weights = 1 / var_comp
  , data = comm_comp_summary_p2)
  
}

#######
## (2) mgcv smooth model to incorporate all predictors (spatio-temporal model)
## For exploration see "spatio_temporal_model.R"
#######

source("GAM_setup.R")

## Subset data for well observed counties
comm_comp_summary_p2_well_observed_m_dh <- comm_comp_summary_p2_m_dh %>% filter(num_lists >= 80)

source("GAM_run.R")

source("GAM_cross_validation.R")

#######
## (3) Calculating species importance (for now without error at the level of the community)
## Skip if running for modeling the effects of uncertainty for the benefit of time
#######

source("check_species_importance.R")

#######
### (4) Gathering results for comparing uncertainty runs
#######

summary_stats <- list(

  med_comp_f    = comm_comp_summary_p2[["med_comp"]]
, med_comp_r    = comm_comp_summary_p2_well_observed[["med_comp"]]
, var_comp_f    = comm_comp_summary_p2[["var_comp"]]
, var_comp_r    = comm_comp_summary_p2_well_observed[["var_comp"]]
, cv_comp_f     = comm_comp_summary_p2[["cv_comp"]]
, cv_comp_r     = comm_comp_summary_p2_well_observed[["cv_comp"]]
, spatio_r.sq_f = summary(comm_comp_spatio_temporal_f[[2]])[["r.sq"]]
, spatio_r.sq_r = summary(comm_comp_spatio_temporal_r[[2]])[["r.sq"]]

)

#######
## (5) Compile
#######
final_stats.out <- list(
  dilut_eff        = list(
    dilut_eff_lm_red
  , dilut_eff_pred_lm_red
    )
, spatio_temporal  = list(
    comm_comp_spatio_temporal_f
  , comm_comp_spatio_temporal_r
  )  
, species_specific = list(
    species_importance_p_f
  , species_importance_p_r
  , true_amplifiers_f
  , true_amplifiers_r
  , true_diluters_f
  , true_diluters_r
  , species_importance_most_observed_f
  , species_importance_most_observed_r
  , species_importance_most_abundant_f
  , species_importance_most_abundant_r
    )
, summary_stats    = summary_stats
  )

########
## (6) Placeholder space for extra stuff
## When running scripts in an automated way this stuff isn't needed, but can be run for 
## full analyses
########

## skip over this stuff if modeling uncertainty
if (test_uncertainty_effect == FALSE) {

how_many_gtg <- data.frame(
  species   = numeric(length(red_all_spec))
, num_lists = numeric(length(red_all_spec))
)

how_many_gtg <- species_to_est_data_with_meta_red %>%
  group_by(Scientific_Name) %>%
  summarize(n())

quantile(unlist(how_many_gtg[, 2]), c(0.05, 0.5, 0.95))

as.data.frame(how_many_gtg)[order(as.data.frame(how_many_gtg[, 2])), ]

how_many_gtg[how_many_gtg$Scientific_Name == "Gymnorhinus_cyanocephalus", ]

}
## only set to true if all spatio tempral aggregations are desired. Otherwise run individual manually
if (run_spat_temp_aggregated == TRUE) {
  
####
## Full model
####
comm_comp_summary  <- comp_prop_sum_R0(
  host_comp_summary  = host_comp_summary
, bird_prop_dat      = samp_data_com
, comm_detect_est    = comm_detect_est
, uncertainty_list   = uncertainty_list
, mosquito_survival  = new_temp_data
, m_to_b_trans_samps = samps_mos_bird
, county_temp_data   = county_temp_data)

## clean
source("community_competence_just_clean.R")

comm_comp_summary_p2_well_observed_full_nu <- comm_comp_summary_p2_well_observed
comm_comp_summary_p2_full_nu               <- comm_comp_summary_p2
#saveRDS(comm_comp_summary_p2_well_observed_full_nu, "comm_comp_summary_p2_well_observed_full_nu.Rds")
#saveRDS(comm_comp_summary_p2_full_nu, file = "comm_comp_summary_p2_full_nu.Rds")

########################
## Mean temperature over space
########################
county_temp_data_temp <- county_temp_data %>%
  group_by(month) %>%
  summarize(temp = mean(temp))
county_temp_data_temp <- left_join(county_temp_data, county_temp_data_temp, by = "month")
county_temp_data_temp <- county_temp_data_temp[, -4]
names(county_temp_data_temp)[4] <- "temp"
county_temp_data_temp <- transform(county_temp_data_temp, temp = round(temp))

comm_comp_summary  <- comp_prop_sum_R0(
  host_comp_summary  = host_comp_summary
, bird_prop_dat      = samp_data_com
, comm_detect_est    = comm_detect_est
, uncertainty_list   = uncertainty_list
, mosquito_survival  = new_temp_data
, m_to_b_trans_samps = samps_mos_bird
, county_temp_data   = county_temp_data_temp)

## clean
source("community_competence_just_clean.R")

comm_comp_summary_p2_well_observed_ts_nu <- comm_comp_summary_p2_well_observed
comm_comp_summary_p2_ts_nu <- comm_comp_summary_p2
#saveRDS(comm_comp_summary_p2_well_observed_ts_nu, "comm_comp_summary_p2_well_observed_ts_nu.Rds")
#saveRDS(comm_comp_summary_p2_ts_nu, file = "comm_comp_summary_p2_ts_nu.Rds")

########################
## Mean temperature over time
########################
county_temp_data_temp <- county_temp_data[county_temp_data$month == "5" | county_temp_data$month == "10", ] %>%
  group_by(county) %>%
  summarize(temp = mean(temp))
county_temp_data_temp <- left_join(county_temp_data, county_temp_data_temp, by = "county")
county_temp_data_temp <- county_temp_data_temp[, -4]
names(county_temp_data_temp)[4] <- "temp"
county_temp_data_temp <- transform(county_temp_data_temp, temp = round(temp))

comm_comp_summary  <- comp_prop_sum_R0(
  host_comp_summary  = host_comp_summary
, bird_prop_dat      = samp_data_com
, comm_detect_est    = comm_detect_est
, uncertainty_list   = uncertainty_list
, mosquito_survival  = new_temp_data
, m_to_b_trans_samps = samps_mos_bird
, county_temp_data   = county_temp_data_temp)

## clean
source("community_competence_just_clean.R")

comm_comp_summary_p2_well_observed_tt_nu <- comm_comp_summary_p2_well_observed
comm_comp_summary_p2_tt_nu <- comm_comp_summary_p2
#saveRDS(comm_comp_summary_p2_well_observed_tt_nu, "comm_comp_summary_p2_well_observed_tt_nu.Rds")
#saveRDS(comm_comp_summary_p2_tt_nu, file = "comm_comp_summary_p2_tt_nu.Rds")

########################
## Mean bird community over space
########################
samp_data_com_temp <- samp_data_com %>%
  group_by(month, phylo_name) %>%
  summarize(ebird_prior = mean(ebird_prior))
## create an empty data frame with each county x month x year combination, but with empty 
 ## scientific name slot. The problem is that this requires an insane amount of RAM
temp_ebird_samples <- samp_data_com %>% dplyr::select(county, month, year)
temp_ebird_samples <- temp_ebird_samples[-which(duplicated(temp_ebird_samples) == TRUE), ]
## To save RAM, just run the model with the well sampled dataset
temp_ebird_samples <- comm_comp_summary_p2_well_observed_full %>% dplyr::select(county, month, year)

## sort samp data com temp by rarity of birds
samp_data_com_temp <- left_join(samp_data_com_temp, temp_ebird_samples, by = c("month"))

comm_comp_summary  <- comp_prop_sum_R0(
  host_comp_summary  = host_comp_summary
, bird_prop_dat      = samp_data_com_temp
, comm_detect_est    = comm_detect_est
, uncertainty_list   = uncertainty_list
, mosquito_survival  = new_temp_data
, m_to_b_trans_samps = samps_mos_bird
, county_temp_data   = county_temp_data)

## clean
source("community_competence_just_clean.R")

comm_comp_summary_p2_well_observed_bcs <- comm_comp_summary_p2_well_observed
comm_comp_summary_p2_bcs <- comm_comp_summary_p2
#saveRDS(comm_comp_summary_p2_well_observed_bcs, "comm_comp_summary_p2_well_observed_bcs.Rds")
#saveRDS(comm_comp_summary_p2_bcs, "comm_comp_summary_p2_bcs.Rds")

########################
## Mean bird community over time
########################
samp_data_com_temp <- samp_data_com %>%
  group_by(county, phylo_name) %>%
  summarize(ebird_prior = mean(ebird_prior))

samp_data_com_temp <- samp_data_com_temp %>% filter(county != "")

## To save RAM, just run the model with the well sampled dataset
temp_ebird_samples <- comm_comp_summary_p2_well_observed_full %>% dplyr::select(county, month, year)
## sort samp data com temp by rarity of birds
samp_data_com_temp <- left_join(samp_data_com_temp, temp_ebird_samples, by = c("county"))
## remove the rows where month or year is NA
samp_data_com_temp <- samp_data_com_temp %>% filter(!is.na(month)) %>% filter(!is.na(year))

comm_comp_summary  <- comp_prop_sum_R0(
  host_comp_summary  = host_comp_summary
, bird_prop_dat      = samp_data_com_temp
, comm_detect_est    = comm_detect_est
, uncertainty_list   = uncertainty_list
, mosquito_survival  = new_temp_data
, m_to_b_trans_samps = samps_mos_bird
, county_temp_data   = county_temp_data)

## clean
source("community_competence_just_clean.R")

## Summarize from here with the code in summary_stats_cleaned.R 
comm_comp_summary_p2_well_observed_bct <- comm_comp_summary_p2_well_observed
comm_comp_summary_p2_bct <- comm_comp_summary_p2
#saveRDS(comm_comp_summary_p2_well_observed_bct, "comm_comp_summary_p2_well_observed_bct.Rds")
#saveRDS(comm_comp_summary_p2_bct, "comm_comp_summary_p2_bct.Rds")

################################################
## Mean bird community and mean temperature over time and space
################################################
samp_data_com_temp <- samp_data_com %>%
  group_by(phylo_name) %>%
  summarize(ebird_prior = mean(ebird_prior))

## To save RAM, just run the model with the well sampled dataset
temp_ebird_samples <- comm_comp_summary_p2_well_observed %>% dplyr::select(county, month, year)
temp_ebird_samples <- temp_ebird_samples[rep(seq(1:nrow(temp_ebird_samples)), each = nrow(samp_data_com_temp)), ]

## Add in the birds to each community. Not great to use 
samp_data_com_temp <- transform(temp_ebird_samples
  , phylo_name  = rep(samp_data_com_temp$phylo_name, nrow(comm_comp_summary_p2_well_observed))
  , ebird_prior = rep(samp_data_com_temp$ebird_prior, nrow(comm_comp_summary_p2_well_observed)))

samp_data_com_temp <- samp_data_com_temp[1:nrow(comm_bite_est), ]

## Can do a trick here, where because I only have one county and month that is identical with all of the other
## counties and months, I can just run the model with a single community and then duplicate the data 
## a number of times so that the number of communities in each model are identical for comparison of the 
## MSE

## Now the temp, averaged over space and time
county_temp_data_temp <- county_temp_data %>% summarize(temp = mean(temp))
county_temp_data_temp <- transform(county_temp_data, temp = county_temp_data_temp$temp)
county_temp_data_temp <- transform(county_temp_data_temp, temp = round(temp))

comm_comp_summary  <- comp_prop_sum_R0(
  host_comp_summary  = host_comp_summary
, bird_prop_dat      = samp_data_com_temp
, comm_detect_est    = comm_detect_est
, uncertainty_list   = uncertainty_list
, mosquito_survival  = new_temp_data
, m_to_b_trans_samps = samps_mos_bird
, county_temp_data   = county_temp_data_temp
, print_prog         = FALSE) ## Runs in < 0.3 second per, so a waste to print

## last step is to just replicate the med_comp value predicted for this average temp
 ## average average bird community community equal to the number of well sampled communities
comm_comp_summary_p2_well_observed$med_comp <- comm_comp_summary_p2_well_observed$med_comp[1]

## clean
source("community_competence_just_clean.R")

## Summarize from here with the code in summary_stats_cleaned.R 
comm_comp_summary_p2_well_observed_mean_model <- comm_comp_summary_p2_well_observed
comm_comp_summary_p2_mean_model <- comm_comp_summary_p2
#saveRDS(comm_comp_summary_p2_well_observed_mean_model, "comm_comp_summary_p2_well_observed_mean_model.Rds")
#saveRDS(comm_comp_summary_p2_mean_model, "comm_comp_summary_p2_mean_model.Rds")
  
}
