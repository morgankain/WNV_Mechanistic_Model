################################
### Check species importance ###
################################

if (test_species_importance == TRUE) {

host_comp_summary_m <- melt(host_comp_summary, c("day", "log_dose", "species", "model", "outcome"))

## organizing species competence for plots
host_comp_summary_p1 <- host_comp_summary_m %>% 
  filter(outcome == "realized_transmission") %>%
  group_by(species, variable) %>%
  summarize(tot_trans = sum(value))

host_comp_summary_p2 <- host_comp_summary_m %>% 
  filter(outcome == "biting_pref") %>%
  group_by(species, variable) %>%
  summarize(bite_pref = mean(value))

## Quantiles
host_comp_summary_p1 <- host_comp_summary_p1 %>% 
  group_by(species) %>%
  summarize(med_comp = quantile(tot_trans, probs = 0.5))

host_comp_summary_p2 <- host_comp_summary_p2 %>% 
  group_by(species) %>%
  summarize(bite_pref = quantile(bite_pref, probs = 0.5))

host_comp_summary_p1 <- merge(host_comp_summary_p1, host_comp_summary_p2)

comm_detect_est_m <- melt(comm_detect_est, c("species", "body_size", "model"))

comm_detect_est_p1 <- comm_detect_est_m %>% 
  group_by(species) %>%
  summarize(detectability = quantile(value, probs = 0.5))

host_comp_summary_p1 <- merge(host_comp_summary_p1, comm_detect_est_p1)

  ## merge effort data with community data
species_to_est_data_with_meta <- merge(
    samp_data_com
  , effort_metadata_com[, c(1, 2, 3, 6)]
  , by = c("county", "month", "year"), all = TRUE)

  ## number of communities with species X
species_to_est_data_with_meta_red <- species_to_est_data_with_meta %>%
  filter(num_lists >= 80)

species_to_est_data_with_meta_red <- droplevels(species_to_est_data_with_meta_red)
red_all_spec <- unique(species_to_est_data_with_meta_red[["phylo_name"]])

## see bottom of this script for creating species_to_est_data_with_meta_red
 ## For full data set use:     species_to_est_data_with_meta
 ## For reduced data set use:  species_to_est_data_with_meta_red
species_importance <- check_species_importance_no_uncer_R0(
   host_comp_summary  = host_comp_summary_p1
 , bird_prop_dat      = species_to_est_data_with_meta_red
 , comm_detect_est    = host_comp_summary_p1[, c(1, 4)]
 , nsamps             = ifelse(no_uncer == TRUE, 1, uncertainty_list[["samps"]])
 , mosquito_survival  = mosq_surv
 , m_to_b_trans_samps = samps_mos_bird)

rm(host_comp_summary_m); gc()

## Add meta-data to species importance: total number of lists including each species
species_importance <- species_importance %>%
  group_by(phylo_name) %>%
  mutate(n_obs = n())

## Add meta-data to species importance: total number of lists submitted in each County:Month:Year
species_importance <- merge(
    species_importance
  , effort_metadata_com[, c(1, 2, 3, 6)]
  , by = c("county", "month", "year"))

## Add meta-data to species importance: total number of observations of each bird
total_species_obs <- samp_data_com %>%
  group_by(phylo_name) %>%
  summarize(total_obs = sum(ebird_prior))

species_importance <-  merge(
    species_importance
  , total_species_obs
  , by = c("phylo_name"), all = TRUE)

## Subset the data to well sampled counties, months, and years 
  species_importance_well_observed <- species_importance %>%
  group_by(phylo_name) %>%
  filter(num_lists >= 80)      ## Total # of lists

species_importance_p_f <- 
  species_importance %>%
  group_by(phylo_name) %>%
  summarize(
    med_dilut_effect        = quantile(spec_dilut_amp_substitutive, probs = 0.5, na.rm = TRUE)
  , min_dilut_effect        = quantile(spec_dilut_amp_substitutive, probs = 0.05, na.rm = TRUE)
  , min_dilut_effect_narrow = quantile(spec_dilut_amp_substitutive, probs = 0.20, na.rm = TRUE)
  , max_dilut_effect        = quantile(spec_dilut_amp_substitutive, probs = 0.95, na.rm = TRUE)
  , max_dilut_effect_narrow = quantile(spec_dilut_amp_substitutive, probs = 0.80, na.rm = TRUE)
  , mean_dilut_effect       = mean(spec_dilut_amp_substitutive, na.rm = TRUE)
  , max_true_dilut_effect   = max(spec_dilut_amp_substitutive, na.rm = TRUE)
  , min_true_dilut_effect   = min(spec_dilut_amp_substitutive, na.rm = TRUE)
  , sum_dilut_effect        = sum(spec_dilut_amp_substitutive, na.rm = TRUE)
  , range_dilut_effect      = (max(spec_dilut_amp_substitutive) - 1) + (1 - min(spec_dilut_amp_substitutive))
  , med_rank                = quantile(spec_rank, probs = 0.5, na.rm = TRUE)
  , min_rank                = quantile(spec_rank, probs = 0.05, na.rm = TRUE)
  , max_rank                = quantile(spec_rank, probs = 0.95, na.rm = TRUE)
  , num_lists               = mean(num_lists)
  , n_obs                   = mean(n_obs)
  , total_obs               = mean(total_obs)
  , n_com                   = n()
  , n_com_above_1           = length(which(full_comm_R0 > 1))
  , n_com_below_1           = length(which(full_comm_R0 < 1))
  , over_1_below_1          = sum(over_1_below_1)
  , below_1_over_1          = sum(below_1_over_1)
  , drop_from_above_1       = sum(over_1_below_1) / length(which(full_comm_R0 > 1))
  , move_from_below_1       = sum(below_1_over_1) / length(which(full_comm_R0 < 1)))

species_importance_p_r <- 
  species_importance_well_observed %>%
  group_by(phylo_name) %>%
  summarize(
    med_dilut_effect        = quantile(spec_dilut_amp_substitutive, probs = 0.5, na.rm = TRUE)
  , min_dilut_effect        = quantile(spec_dilut_amp_substitutive, probs = 0.05, na.rm = TRUE)
  , min_dilut_effect_narrow = quantile(spec_dilut_amp_substitutive, probs = 0.20, na.rm = TRUE)
  , max_dilut_effect        = quantile(spec_dilut_amp_substitutive, probs = 0.95, na.rm = TRUE)
  , max_dilut_effect_narrow = quantile(spec_dilut_amp_substitutive, probs = 0.80, na.rm = TRUE)
  , mean_dilut_effect       = mean(spec_dilut_amp_substitutive, na.rm = TRUE)
  , max_true_dilut_effect   = max(spec_dilut_amp_substitutive, na.rm = TRUE)
  , min_true_dilut_effect   = min(spec_dilut_amp_substitutive, na.rm = TRUE)
  , sum_dilut_effect        = sum(spec_dilut_amp_substitutive, na.rm = TRUE)
  , range_dilut_effect      = (max(spec_dilut_amp_substitutive) - 1) + (1 - min(spec_dilut_amp_substitutive))
  , med_rank                = quantile(spec_rank, probs = 0.5, na.rm = TRUE)
  , min_rank                = quantile(spec_rank, probs = 0.05, na.rm = TRUE)
  , max_rank                = quantile(spec_rank, probs = 0.95, na.rm = TRUE)
  , num_lists               = mean(num_lists)
  , n_obs                   = mean(n_obs)
  , total_obs               = mean(total_obs)
  , n_com                   = n()
  , n_com_above_1           = length(which(full_comm_R0 > 1))
  , n_com_below_1           = length(which(full_comm_R0 < 1))
  , over_1_below_1          = sum(over_1_below_1)
  , below_1_over_1          = sum(below_1_over_1)
  , drop_from_above_1       = sum(over_1_below_1) / length(which(full_comm_R0 > 1))
  , move_from_below_1       = sum(below_1_over_1) / length(which(full_comm_R0 < 1)))

#######
## Organize and summarize species importance
#######

######
## Species that have the highest median for an amplification effect
######
true_amplifiers_f <- species_importance_p_f[order(species_importance_p_f[["max_true_dilut_effect"]], decreasing = FALSE), ]
true_amplifiers_f <- true_amplifiers_f[order(true_amplifiers_f[["med_dilut_effect"]], decreasing = FALSE), ]
true_amplifiers_f[["phylo_name"]] <- factor(true_amplifiers_f[["phylo_name"]]
  , levels = rev(true_amplifiers_f[["phylo_name"]]))
true_amplifiers_f <- true_amplifiers_f[1:15, ]

true_amplifiers_r <- species_importance_p_r[order(species_importance_p_r[["max_true_dilut_effect"]], decreasing = FALSE), ]
true_amplifiers_r <- true_amplifiers_r[order(true_amplifiers_r[["med_dilut_effect"]], decreasing = FALSE), ]
true_amplifiers_r[["phylo_name"]] <- factor(true_amplifiers_r[["phylo_name"]]
  , levels = rev(true_amplifiers_r[["phylo_name"]]))
true_amplifiers_r <- true_amplifiers_r[1:15, ]

######
## Species that have the highest median for a dilution effect
######
true_diluters_f <- species_importance_p_f[order(species_importance_p_f[["min_true_dilut_effect"]], decreasing = TRUE), ]
true_diluters_f <- true_diluters_f[order(true_diluters_f[["med_dilut_effect"]], decreasing = TRUE), ]
true_diluters_f[["phylo_name"]] <- factor(true_diluters_f[["phylo_name"]]
  , levels = rev(true_diluters_f[["phylo_name"]]))
true_diluters_f <- true_diluters_f[1:15, ]

true_diluters_r <- species_importance_p_r[order(species_importance_p_r[["min_true_dilut_effect"]], decreasing = TRUE), ]
true_diluters_r <- true_diluters_r[order(true_diluters_r[["med_dilut_effect"]], decreasing = TRUE), ]
true_diluters_r[["phylo_name"]] <- factor(true_diluters_r[["phylo_name"]]
  , levels = rev(true_diluters_r[["phylo_name"]]))
true_diluters_r <- true_diluters_r[1:15, ]

######
## Most observed species across Texas
######
species_importance_p_f <- species_importance_p_f[order(species_importance_p_f[["n_obs"]], decreasing = TRUE), ]
species_importance_p_f[["phylo_name"]] <- factor(species_importance_p_f[["phylo_name"]]
  , levels = rev(species_importance_p_f[["phylo_name"]]))
species_importance_most_observed_f <- species_importance_p_f[1:15, ]

species_importance_p_r <- species_importance_p_r[order(species_importance_p_r[["n_obs"]], decreasing = TRUE), ]
species_importance_p_r[["phylo_name"]] <- factor(species_importance_p_r[["phylo_name"]]
  , levels = rev(species_importance_p_r[["phylo_name"]]))
species_importance_most_observed_r <- species_importance_p_r[1:15, ]

######
## Most abundant species across Texas
######
species_importance_p_f <- species_importance_p_f[order(species_importance_p_f[["total_obs"]], decreasing = TRUE), ]
species_importance_p_f[["phylo_name"]] <- factor(species_importance_p_f[["phylo_name"]]
  , levels = rev(species_importance_p_f[["phylo_name"]]))
species_importance_most_abundant_f <- species_importance_p_f[1:15, ]

species_importance_p_r <- species_importance_p_r[order(species_importance_p_r[["total_obs"]], decreasing = TRUE), ]
species_importance_p_r[["phylo_name"]] <- factor(species_importance_p_r[["phylo_name"]]
  , levels = rev(species_importance_p_r[["phylo_name"]]))
species_importance_most_abundant_r <- species_importance_p_r[1:15, ]

} else {
  ## if modeling the effects of uncertainty skip these
  species_importance_p_f              <- NULL
  species_importance_p_r              <- NULL
  true_amplifiers_f                   <- NULL
  true_amplifiers_r                   <- NULL
  true_diluters_f                     <- NULL
  true_diluters_r                     <- NULL
  species_importance_most_observed_f  <- NULL
  species_importance_most_observed_r  <- NULL
  species_importance_most_abundant_f  <- NULL
  species_importance_most_abundant_r  <- NULL
  
}
