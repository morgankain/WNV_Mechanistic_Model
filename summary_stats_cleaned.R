######################################################################
### Calculate some summary stats, cleaned for automated script run ###
######################################################################

#######
## First, test some assumptions of the dilution effect
## Second, fit a spatio-temporal model
## Third, calculate the importance of species in a few different ways
#######

## Subset data for well observed counties
comm_comp_summary_p2_well_observed <- comm_comp_summary_p2 %>%
  filter(num_lists >= 80)

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
  , data.frame(detectability = comm_detect_est[, -c(1:3)])
)

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

## combine adjusted communities with estimated species physioligcal competence
dilut_assump_full <- left_join(host_comp_summary_p2[, c(1, 2, 4)]
  , species_to_est_data_p)

dilut_assump_red <- left_join(host_comp_summary_p2[, c(1, 2, 4)]
  , species_to_est_data_p2)

if (no_uncer == TRUE) {

dilut_eff_lm_red           <- lm(med_comp ~ log(med_obs)
  , data = dilut_assump_red)

dilut_eff_pred_lm_red      <- lm(med_comp ~ log(meannum_spec)
  , data = comm_comp_summary_p2_well_observed)

} else {
  
dilut_eff_lm_red      <- lm(med_comp ~ log(med_obs)
  , weights = 1 / var_comp
  , data = dilut_assump_red)

dilut_eff_pred_lm_red <- lm(med_comp ~ log(meannum_spec)
  , weights = 1 / var_comp
  , data = comm_comp_summary_p2_well_observed)

dilut_eff_pred_lm <- lm(med_comp ~ log(meannum_spec)
  , weights = 1 / var_comp
  , data = comm_comp_summary_p2)
  
}

gg_diult_eff <- ggplot(dilut_assump_full, aes(med_obs, med_comp)) + 
  geom_point(lwd = 1) + 
  geom_smooth(lwd = 1, method = "lm") +
  scale_x_log10(breaks = c(1e-05, 1e-04, 1e-03, 1e-02)) +
  xlab("Median Species Abundance") + 
  ylab("Species Physiological Competence")

#######
## (2) mgcv smooth model to incorporate all predictors (spatio-temporal model)
## For exploration see "spatio_temporal_model.R"
#######

source("GAM_setup.R")

## Subset data for well observed counties
comm_comp_summary_p2_well_observed_m_dh <- comm_comp_summary_p2_m_dh %>%
  filter(num_lists >= 80)

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
## skip over this stuff if modeling uncertainty
## When running scripts in an automated way this stuff isn't needed
########

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

