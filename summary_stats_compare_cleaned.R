########################################################################################
### Automated script for comparing results uncertainty ``on'' vs uncertainty ``off'' ###
########################################################################################

## Goal here is to loop over all of the files in the folder, ignoring the ``master'' file of all uncerainty on
  ## and compare the metrics in each loaded file with the ``master'' file

# med_comp_f    
# med_comp_r  
# var_comp_f   
# var_comp_r   
# cv_comp_f  
# cv_comp_r     
# spatio_r.sq_f 
# spatio_r.sq_r 
# dilut_CI   

## find all of the saved file names in the model_out folder
  ## bit of an ugly way to igonre the folder with the master list
## ignore the subfolder and all uncertainty off for now 
summary_stats_file_list   <- list.files("model_out/")[-c(1, 2)] 

uncer_results <- data.frame(
    model = summary_stats_file_list
,   metric = rep(c(
  "med_comp_f"    
, "med_comp_r"  
, "var_comp_f"   
, "var_comp_r"   
, "cv_comp_f" 
, "cv_comp_r"     
, "spatio_r.sq_f" 
, "spatio_r.sq_r" 
, "dilut_CI" )
  , each = length(summary_stats_file_list))
  , result = numeric(length(summary_stats_file_list) * 9))

all_uncer_model <- readRDS("model_out/_all_uncer/Texas_consensus_1000_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_1_0_1")
no_uncer_model <- readRDS("model_out/_no_uncer/Texas_consensus_1000_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_1")

for (i in seq_along(summary_stats_file_list)) {
  
  temp_model_compare <- summary_stats_file_list[i]
  
  one_off_uncer      <- readRDS(paste("model_out/", temp_model_compare, sep = ""))
  
  ## all metrics for model to compare for this run
  results_compare <- data.frame(
    model = temp_model_compare
,   metric = c(
  "med_comp_f"    
, "med_comp_r"  
, "var_comp_f"   
, "var_comp_r"   
, "cv_comp_f" 
, "cv_comp_r"     
, "spatio_r.sq_f" 
, "spatio_r.sq_r" 
, "dilut_CI" 
  )
, result = c(
  median(all_uncer_model[[3]]$summary_stats[["med_comp_f"]] / one_off_uncer[[3]]$summary_stats[["med_comp_f"]])
, median(all_uncer_model[[3]]$summary_stats[["med_comp_r"]] / one_off_uncer[[3]]$summary_stats[["med_comp_r"]])
, median(all_uncer_model[[3]]$summary_stats[["var_comp_f"]] / one_off_uncer[[3]]$summary_stats[["var_comp_f"]])
, median(all_uncer_model[[3]]$summary_stats[["var_comp_r"]] / one_off_uncer[[3]]$summary_stats[["var_comp_r"]])
, median(all_uncer_model[[3]]$summary_stats[["cv_comp_f"]]  / one_off_uncer[[3]]$summary_stats[["cv_comp_f"]])
, median(all_uncer_model[[3]]$summary_stats[["cv_comp_r"]]  / one_off_uncer[[3]]$summary_stats[["cv_comp_r"]])
, (all_uncer_model[[3]]$summary_stats[["spatio_r.sq_f"]]    / one_off_uncer[[3]]$summary_stats[["spatio_r.sq_f"]])
, (all_uncer_model[[3]]$summary_stats[["spatio_r.sq_r"]]    / one_off_uncer[[3]]$summary_stats[["spatio_r.sq_r"]])
, (all_uncer_model[[3]]$summary_stats[["dilut_CI"]]         / one_off_uncer[[3]]$summary_stats[["dilut_CI"]])
)
)
 
  uncer_results[uncer_results$model == temp_model_compare, ] <- results_compare

  print(i)
}

uncer_results[uncer_results$model == "Texas_consensus_1000_1_1_1_1_1_1_0_1_1_1_1_1_1_1_1_1_0_1", ]

uncer_results[uncer_results$model == "Texas_consensus_1000_1_1_1_1_1_1_0_1_1_1_1_1_1_1_1_1_0_1", ][c(4, 6, 8, 9), ]
