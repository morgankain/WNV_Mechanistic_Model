###################################################################################################
### finding the tips numbers from the full phylogeny that correspond to the species that I need ###
###################################################################################################

specs_to_est_titer_surv <- match_scinames(
  needed_species = spec_est_vec
, have_species   = needed_species
, phylo          = sing_tree
)

specs_to_est_bite <- match_scinames(
  needed_species = spec_est_vec
, have_species   = data.frame(Scientific.Name = bite_pref_res[ , 3])   ## switch the order of the columns for the function
, phylo          = sing_tree
)