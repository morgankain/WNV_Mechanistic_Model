#########################################
### Stan model for biting preferences ###
#########################################

if (load_stan_bite_pref_res == TRUE & file.exists("saved_output/mosq_bite.Rds")) {
 
mosq_bite_rds <- readRDS("saved_output/mosq_bite.Rds")

bite_pref_res        <- mosq_bite_rds[[1]]
mosq_bite.data       <- mosq_bite_rds[[2]]
mosq_bite_model_out  <- mosq_bite_rds[[3]]
mosq_bite.summary    <- mosq_bite_rds[[4]]
rm(mosq_bite_rds)
 
names(bite_pref_res)[1:3] <- c("Bites", "Preference", "Scientific_Name")

  } else {

  ## Data
mosq_bite.data <- list(
       "N"            = length(spec_prop_comp[["Proportion"]])
    ,  "hosts_prop"   = host_species
    ,  "prop"         = spec_prop_comp[["Proportion"]]
    ,  "bites"        = spec_bites[["Bites"]]
    ,  "flat_alpha_p" = rep(1, length(spec_prop_comp[["Proportion"]]))
    ,  "prev_obs_p"   = prev_obs_p[["Ebird_Prior"]]
    ,  "prev_obs_b"   = prev_obs_b)

  ## Run Model
mosq_bite_model_out <- stan(
  file    = "stan/mosq_bite.stan"
, data    = mosq_bite.data
, iter    = 10000
, thin    = 2
, warmup  = 1000
, refresh = max(10000/100, 1)
, control = list(max_treedepth = 12, adapt_delta = .92)
, chains  = 4)

## Pleasant way to look at convergence of the model
# launch_shinystan(mosq_bite_model_out)

## organize model output
detach("package:tidyr", unload = TRUE)
samps_mosq_bite          <- extract(mosq_bite_model_out, permuted = FALSE)
library(tidyr)
tidy_mosq_bite           <- tidy(mosq_bite_model_out)
mosq_bite_model_out_summ <- summary(mosq_bite_model_out)
mosq_bite                <- mosq_bite_model_out_summ$summary

## summarize data and estimates
mosq_bite.summary <- data.frame(
  data             = c(mosq_bite.data[["prop"]], mosq_bite.data$bites)
, estimate         = c(
  mosq_bite[grep("theta_p", dimnames(mosq_bite)[[1]])
    , grep("50%", dimnames(mosq_bite)[[2]])]
, mosq_bite[grep("bite_pref", dimnames(mosq_bite)[[1]])
    , grep("50%", dimnames(mosq_bite)[[2]])])
, outcome          = c(rep("bird_prop", 150), rep("bite_preference", 150))
, species          = rep(mosq_bite_prop_trial[["Scientific_Name"]], 2)
, prop_prior       = mosq_bite.data[["flat_alpha_p"]]
, prop_prior_count = mosq_bite.data[["prev_obs_p"]]
, prop_prior_bite  = prev_obs_b)

## Sets up the data frame for the phylo model for biting
bite_pref_res        <- mosq_bite.summary %>% dplyr::select(data, estimate, species) 
bite_pref_res        <- bite_pref_res[grep("bite_pref", dimnames(mosq_bite)[[1]]), ]

### Adjust pref for Poisson lme4 model: multiplying by 1000 and rounding will have minimal effect
bite_pref_res  <- transform(bite_pref_res, Scaled_Preference = round(estimate * 1000))

## Add columns for uncertainty (# draws = uncertainty_list[["samps"]])
 ## take 1/4 of samples from each chain
  ## Need biting pref samples
samps_mosq_bite <- samps_mosq_bite[,, grep("bite_pref", dimnames(mosq_bite)[[1]])]

for (i in 1:4) {
  
  if (i == 1) {
    bites_pred_samp_mat <- t(sample_n(as.data.frame(samps_mosq_bite[,i,])
    , round(uncertainty_list[["samps"]] / 4))) * 1000
  } else {
    bites_pred_samp_mat <- cbind(bites_pred_samp_mat
      , t(sample_n(as.data.frame(samps_mosq_bite[,i,])
    , round(uncertainty_list[["samps"]] / 4))) * 1000)
  }
}

## Combine median with all samples
bite_pref_res <- cbind(bite_pref_res, bites_pred_samp_mat)

  saveRDS(
    list(
      bite_pref_res
    , mosq_bite.data
    , mosq_bite_model_out
    , mosq_bite.summary
    )
    , file = "saved_output/mosq_bite.Rds")
  
  names(bite_pref_res)[1:3] <- c("Bites", "Preference", "Scientific_Name")
  
  }
