#############################################################
### Predict responses for species outside of the data set ###
### Predict mosquito biting preferences for new species   ###
#############################################################

## If previous responses are to be loaded, and previous responses exist for the combination of uncertainty
  ## being considered for the given run, load the previously saved results instead of going through the whole function
## Also, double check if all species are present. If not, default back to running the function
if (used_saved_responses_bite == TRUE) {
  
  ## determine the file name for the current uncertainty run
  comm_bite_est <- response_load_logcial(which_model = "bite", uncertainty_list)
  
  ## If there are no saved responses for the current combination of uncertainty, run the function
  if(is.null(comm_bite_est)) {

comm_bite_est <- outside_spec_est(
  phylo                  = MyTree
, spec_with_dat          = data.frame(bite_pref_res[["Scientific_Name"]]) 
, have_tips              = specs_to_est_bite[["have_tips"]]
, need_tips              = specs_to_est_bite[["needed_tips"]]
, missing_spec           = specs_to_est_bite[["missing_spec"]]
, which_species_no_phylo = which_species_no_phylo
, data.full              = bite_pref_res
, community_spec_bs      = community_spec_bs
, which_model            = "bite"
, outside                = TRUE
, use_saved_responses    = used_saved_responses_bite
, write_resp             = TRUE
, uncertainty_list       = uncertainty_list
, other_responses        = NULL)
}

  } else {
  
comm_bite_est <- outside_spec_est(
  phylo                  = MyTree
, spec_with_dat          = data.frame(bite_pref_res[["Scientific_Name"]]) 
, have_tips              = specs_to_est_bite[["have_tips"]]
, need_tips              = specs_to_est_bite[["needed_tips"]]
, missing_spec           = specs_to_est_bite[["missing_spec"]]
, which_species_no_phylo = which_species_no_phylo
, data.full              = bite_pref_res
, community_spec_bs      = community_spec_bs
, which_model            = "bite"
, outside                = TRUE
, use_saved_responses    = used_saved_responses_bite
, write_resp             = TRUE
, uncertainty_list       = uncertainty_list
, other_responses        = NULL)
     
}

## Remove a few species that were duplicated (having two ebird names but only a single phylo name,
  ## that is they were combined in the most recent phylogeny)
comm_bite_est <- comm_bite_est[!duplicated(comm_bite_est[["species"]]), ]

## remove those species that don't show up in the data
comm_bite_est <- comm_bite_est[!is.na(match(comm_bite_est[["species"]]
  , samp_data_com[["phylo_name"]])), ]

## If there is no uncertainty at all, remove all but the first column
if (no_uncer == TRUE) {
  comm_bite_est <- comm_bite_est[, 1:3]
}

