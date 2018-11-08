#######################################
### Control what is written to disk ###
#######################################

  ## Build a name to save the results
model_out.filename <- paste(
  which_state
, which_phy
, paste(unlist(uncertainty_list), collapse = "_")
, sep = "_")

saveRDS(list(
  final_spec_comp.out = host_comp_summary
, final_com_comp.out  = comm_comp_summary
, final_stats.out     = final_stats.out)
, paste("model_out"
  , paste(model_out.filename, which_uncertainty_run, sep = "__")
  , sep = "/"))
