#############################################################
### Predict responses for species outside of the data set ###
### Predict titer and survival for new species            ###
#############################################################

## Idea here is that I produce a phylogeny that has the 47 species for which I have data + a single
  ## bird from the cook county data set for which I need to predict, giving a full phylogeny of 48 species

## If previous responses are to be loaded, and previous responses exist for the combination of uncertainty
  ## being considered for the given run, load the previously saved results instead of going through the whole function
## Also, double check if all species are present. If not, default back to running the function
if (used_saved_responses_titer == TRUE) {
  
  ## detemine the file name for the current uncertainty run
  comm_titer_est    <- response_load_logcial(which_model = "titer", uncertainty_list)
  comm_survival_est <- response_load_logcial(which_model = "survival", uncertainty_list)
  
  ## If there are no saved responses for the current combination of uncertainty,
    ## revert to running from scratch 
  ## Survival model requires a fully run (and not subset titer model) so rerun titer if survival is also
    ## not already saved
  if (is.null(comm_titer_est) & is.null(comm_survival_est)) {

comm_titer_est <- outside_spec_est(
  phylo                   = MyTree
 , spec_with_dat          = needed_species
 , have_tips              = specs_to_est_titer_surv[["have_tips"]]
 , need_tips              = specs_to_est_titer_surv[["needed_tips"]]
 , missing_spec           = specs_to_est_titer_surv[["missing_spec"]]
 , which_species_no_phylo = which_species_no_phylo
 , data.full              = titercurves_reduced
 , community_spec_bs      = community_spec_bs
 , which_model            = "titer"
 , outside                = TRUE
 , use_saved_responses    = used_saved_responses_titer
 , write_resp             = TRUE
 , uncertainty_list       = uncertainty_list
 , other_responses        = NULL)
  }
  
  } else {
  
comm_titer_est <- outside_spec_est(
  phylo                   = MyTree
 , spec_with_dat          = needed_species
 , have_tips              = specs_to_est_titer_surv[["have_tips"]]
 , need_tips              = specs_to_est_titer_surv[["needed_tips"]]
 , missing_spec           = specs_to_est_titer_surv[["missing_spec"]]
 , which_species_no_phylo = which_species_no_phylo
 , data.full              = titercurves_reduced
 , community_spec_bs      = community_spec_bs
 , which_model            = "titer"
 , outside                = TRUE
 , use_saved_responses    = used_saved_responses_titer
 , write_resp             = TRUE
 , uncertainty_list       = uncertainty_list
 , other_responses        = NULL)
    
}

## If previous responses are to be loaded, and previous responses exist for the combination of uncertainty
  ## being considered for the given run, load the previously saved results instead of going through the whole function
## Also, double check if all species are present. If not, default back to running the function
if (used_saved_responses_survival == TRUE) {
  
  ## detemine the file name for the current uncertainty run
  comm_survival_est <- response_load_logcial(which_model = "survival", uncertainty_list)
  
  ## If there are no saved responses for the current combination of uncertainty, run the function
  if (is.null(comm_survival_est)) {
    
comm_survival_est <- outside_spec_est(
  phylo                   = MyTree
 , spec_with_dat          = needed_species
 , have_tips              = specs_to_est_titer_surv[["have_tips"]]
 , need_tips              = specs_to_est_titer_surv[["needed_tips"]]
 , missing_spec           = specs_to_est_titer_surv[["missing_spec"]]
 , which_species_no_phylo = which_species_no_phylo
 , data.full              = survival_reduced
 , community_spec_bs      = community_spec_bs
 , which_model            = "survival"
 , outside                = TRUE
 , use_saved_responses    = used_saved_responses_survival
 , write_resp             = TRUE
 , uncertainty_list       = uncertainty_list
 , other_responses        = comm_titer_est)  
  
}
  } else {
    
  ## obtain estimates for survival
comm_survival_est <- outside_spec_est(
  phylo                   = MyTree
 , spec_with_dat          = needed_species
 , have_tips              = specs_to_est_titer_surv[["have_tips"]]
 , need_tips              = specs_to_est_titer_surv[["needed_tips"]]
 , missing_spec           = specs_to_est_titer_surv[["missing_spec"]]
 , which_species_no_phylo = which_species_no_phylo
 , data.full              = survival_reduced
 , community_spec_bs      = community_spec_bs
 , which_model            = "survival"
 , outside                = TRUE
 , use_saved_responses    = used_saved_responses_survival
 , write_resp             = TRUE
 , uncertainty_list       = uncertainty_list
 , other_responses        = comm_titer_est)  
  
  
}

### remove those species that don't show up in the data
comm_titer_est    <- comm_titer_est[!is.na(match(comm_titer_est[["species"]]
  , samp_data_com[["phylo_name"]])), ]
comm_survival_est <- comm_survival_est[!is.na(match(comm_survival_est[["species"]]
  , samp_data_com[["phylo_name"]])), ]

## Remove a few species that were duplicated (having two ebird names but only a single phylo name,
  ## that is they were combined in the most recent phylogeny)
comm_titer_est    <- comm_titer_est[!duplicated(comm_titer_est[, 1:3]), ]
comm_survival_est <- comm_survival_est[!duplicated(comm_survival_est[, 1:3]), ]

## If there is no uncertainty at all, remove all but the first column
if (no_uncer == TRUE) {
  comm_titer_est    <- comm_titer_est[, 1:5]
  comm_survival_est <- comm_survival_est[, 1:6]
}
