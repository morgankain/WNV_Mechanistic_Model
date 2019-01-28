#################################################
### Establish what type of uncertainty to use ###
#################################################

uncertainty_list <- list(
  samps                       = 1000
, titer_model.fixed_uncer     = TRUE
, titer_model.phylo_rand      = TRUE
, titer_model.other_rand      = TRUE
, titer_model.phylo_tip       = TRUE
, survival_model.fixed_uncer  = TRUE
, survival_model.phylo_rand   = TRUE
, survival_model.other_rand   = TRUE
, survival_model.phylo_tip    = TRUE
, bite_model.fixed_uncer      = TRUE
, bite_model.phylo_rand       = TRUE
, bite_model.phylo_tip        = TRUE
, detect_model.fixed_uncer    = TRUE
, detect_model.phylo_rand     = TRUE
, detect_model.phylo_tip      = TRUE
, stan_BtoM_model             = TRUE
, stan_MtoB_model             = TRUE
, stan_bite_model             = TRUE)
