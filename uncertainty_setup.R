#################################################
### Establish what type of uncertainty to use ###
#################################################

uncertainty_list <- list(
  samps                       = 1000
, titer_model.fixed_uncer     = FALSE
, titer_model.phylo_rand      = FALSE
, titer_model.other_rand      = FALSE
, titer_model.phylo_tip       = FALSE
, survival_model.fixed_uncer  = FALSE
, survival_model.phylo_rand   = FALSE
, survival_model.other_rand   = FALSE
, survival_model.phylo_tip    = FALSE
, bite_model.fixed_uncer      = FALSE
, bite_model.phylo_rand       = FALSE
, bite_model.phylo_tip        = FALSE
, detect_model.fixed_uncer    = FALSE
, detect_model.phylo_rand     = FALSE
, detect_model.phylo_tip      = FALSE
, stan_BtoM_model             = FALSE
, stan_MtoB_model             = FALSE
, stan_bite_model             = FALSE)
