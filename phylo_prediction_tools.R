###################################
### Functions for all estimates ###
###################################

##########
### lme4 model fits (low level functions)
##########
lme4_titer_model_fit    <- function (data, phylo, phyloZ, nsp) {

## Titer is a Ricker function of Day, where the slope over time (Ricker of Day) is
 ## affected by both the dose a bird receives and by its body size. The slope over Day
  ## is allowed to vary by Citation, infection experiment (unique_line) and by species
   ## With mode data, these grouping variables could also vary by log_dose
  
phylo_lmm(
	  log(Titer) ~ (log(Day) + Day) * log(body_size_s) + Log_Dose
	  + (1 | Citation)
	  + (1 | unique_line)
	  + (1 + (log(Day) + Day) | Scientific_Name) 
		, data    = data
		, phylonm = c("Scientific_Name")
		, phyloZ  = phyloZ
		, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    , REML = TRUE)
  
}
lme4_survival_model_fit <- function (data, phylo, phyloZ, nsp) {

## Idea here is that survival will depend primarily on titer, and secondarily
 ## on initial dose and day (outside of the effects of titer) as well as body size. 
  ## How birds are effected by titer may vary by citation and infection experiment, 
   ## but definitely by species, which is the focus
phylo_glmm(
    cbind(Died, Alive) ~ 
      Titer + log(body_size_s) + Day
    + (1 | Scientific_Name)
    + (1 | Citation)
    + (1 | unique_line)
    	, data    = data
  	  , family  = binomial(link = "cloglog")
  		, phylonm = c("Scientific_Name")
		  , phylo   = phylo
	  	, phyloZ  = phyloZ
		  , control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))

}
lme4_biting_model_fit   <- function (data, phylo, phyloZ, nsp, uncertainty_list) {

  ## If uncertainty is/is not to be propagated from the stan model for biting preference
  if (uncertainty_list[["stan_bite_model"]] == FALSE) {
    
  phylo_glmm(
	  Scaled_Preference ~ (1 | Scientific_Name) 
		, data    = data
    , family  = poisson 
		, phylonm = c("Scientific_Name")
		, phylo   = phylo
		, phyloZ  = phyloZ
		, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))
    
  } else {
    
data_uncer <- melt(subset(data, select = -c(Scaled_Preference)), c("Bites", "Preference", "Scientific_Name"))
names(data_uncer)[5] <- "Scaled_Preference"
data_uncer[["Scaled_Preference"]] <- round(data_uncer[["Scaled_Preference"]])

  phylo_glmm(
	  Scaled_Preference ~ (1 | Scientific_Name) 
		, data    = data_uncer
    , family  = poisson 
		, phylonm = c("Scientific_Name")
		, phylo   = phylo
		, phyloZ  = phyloZ
		, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))

  }

}
lme4_detect_model_fit   <- function (data, phylo, phyloZ, nsp) {

  phylo_lmm(
	  log(Detection_Distance) ~ log(body_size_s) +
      (1 | Scientific_Name) 
		, data    = data
		, phylonm = c("Scientific_Name")
		, phylo   = phylo
		, phyloZ  = phyloZ
		, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    , REML    = TRUE)
  
}

#######
### Set up random effect vectors (needed within the functions for extracting random effects from models) (low level functions)
#######
lme4_rand_ef_vec <- function (model_coefs, rand_eff_est, rand_ef_id, rand_ef_lengths, names_vec, lme4_fit, is_sd, phylo_data) {

## All of this is pulling out the "b" (the conditional modes of the random effects of
 ## the fitted lme4 model.)
  
## Estimates and names of the random effects
ran_ef_levels      <- getME(lme4_fit, "flist")
which_rand_names   <- apply(rand_ef_id, 2, function(x) paste(x, collapse = "_"))
ran_ef_terms       <- unlist(lapply(lme4_fit@cnms, function (x) length(x))) 
total_ran_efs      <- length(which_rand_names)
temp_ran_names     <- lapply(ran_ef_levels, function(x) rep(unique(x)))

  if (is_sd == FALSE) {
  
## Extract the levels of the random effect
temp_nam_f        <- vector("character")
temp_ran_eff_f    <- vector("character")
ranef_num_f       <- vector("character")

## Build a data frame with the names, levels and estimates. Set up so that all random 
 ## effect structures are _hypothetically_ supported so that 
for (j in seq_along(temp_ran_names)) {
  ranef_num       <- grep(names(ran_ef_levels)[j], which_rand_names)
  if (names(ran_ef_levels)[j] == "Scientific_Name") {
  temp_nam        <- rep(seq(1, length(phylo_data@Dimnames[[2]]), by = 1), each = ran_ef_terms[j])
  temp_ran_eff    <- rep(which_rand_names[ranef_num], length(phylo_data@Dimnames[[2]]))
  ranef_num       <- rep(ranef_num, length(phylo_data@Dimnames[[2]]))
  } else {
  temp_nam        <- as.character(rep(temp_ran_names[[j]], each = ran_ef_terms[j]))  
  temp_ran_eff    <- rep(which_rand_names[ranef_num], length(temp_nam))
  ranef_num       <- rep(ranef_num, length(temp_ran_names[[j]]))
  }
  temp_nam_f      <- c(temp_nam_f, temp_nam)
  temp_ran_eff_f  <- c(temp_ran_eff_f, temp_ran_eff)
  ranef_num_f     <- c(ranef_num_f, ranef_num)
}

all_ranef <- data.frame(
  which_mod      = temp_nam_f
, cond_mod       = getME(lme4_fit, "b")[, 1]
, which_rand     = temp_ran_eff_f
, which_mod_coef = ranef_num_f)
  
for (i in 1:length(rand_ef_lengths)) {
  model_coefs[[i]]       <- filter(all_ranef, which_mod_coef == i) %>% dplyr::select(cond_mod)
  names(model_coefs)[i]  <- which_rand_names[i]
}

full_vcov       <- suppressWarnings(as.data.frame(summary(lme4_fit)[["varcor"]]))
resid_vcov      <- full_vcov[grep("Residual", full_vcov[["grp"]]), ]

model_coefs[[i + 1]]        <- resid_vcov
names(model_coefs)[[i + 1]] <- "Residual_uncertainty"

for (j in seq_along(temp_ran_names)) {
  if (names(ran_ef_levels)[j] == "Scientific_Name") {
  temp_nam <- "phylo_vcov"
  assign(temp_nam, VarCorr(lme4_fit)[[grep("Scientific_Name", names(VarCorr(lme4_fit)))]])
  } else {
  temp_nam    <- paste(names(ran_ef_levels)[j], "vcov", sep = "_")
  assign(temp_nam, full_vcov[grep(names(ran_ef_levels)[j], full_vcov[["grp"]]), ])
  }
model_coefs[[i + 1 + j]]        <- get(temp_nam)  
names(model_coefs)[[i + 1 + j]] <- temp_nam
}

return(model_coefs)

} else {
  
cond_cov_mat   <- lme4:::condVar(lme4_fit)
phylo_coefs_sd <- vector("list", length = total_ran_efs)

## Just need the phylogenetic random effect here
ran_ef_terms <- apply(matrix(unique(rand_ef_lengths)), 1, function(x) length(which(rand_ef_lengths == x)))
names_vec    <- rep(names(lme4_fit@cnms), ran_ef_terms)

## row and col of the first val for each random effect
start_loc <- cumsum(rand_ef_lengths) - rand_ef_lengths + 1

## row and col of the last val for each random effect
if (length(start_loc) > 1) {
  end_loc <- start_loc[-1] - 1
} else {
  end_loc <- rand_ef_lengths
}

which_sci_rand <- grep("Scientific_Name", which_rand_names)

ranef2 <- rand_eff_est[
  min(start_loc[which_sci_rand]):max(end_loc[which_sci_rand])
, min(start_loc[which_sci_rand]):max(end_loc[which_sci_rand])]

ranef2_info <- data.frame(
  ranef_level    = rep(seq(1, length(phylo_data@Dimnames[[2]]), by = 1), each = ran_ef_terms[grep("Scientific_Name", names(temp_ran_names))])
, ranef          = rep(which_rand_names[which_sci_rand], length(phylo_data@Dimnames[[2]]))
, which_mod_coef = rep(which_sci_rand, length(phylo_data@Dimnames[[2]])))

ranef2 <- cbind(ranef2_info, as.matrix(ranef2))

## Sort out phylo random effect
for (i in which_sci_rand) {
  which_subset          <- which(ranef2[, 3] == i)
  model_coefs[[i]]      <- ranef2[, -c(1, 2, 3)][which_subset, which_subset]
  names(model_coefs)[i] <- paste(paste(rand_ef_id[, i], collapse = "_"), "sd", sep = "_")
}

condvar_branch_array <- array(
  data = 0
  , dim = c(length(which_sci_rand)
  , length(which_sci_rand)
  , rand_ef_lengths[which_sci_rand[1]]))

submat_size_start <- seq(1, nrow(ranef2), by = dim(condvar_branch_array)[1])
submat_size_end   <- seq(dim(condvar_branch_array)[1], nrow(ranef2), by = dim(condvar_branch_array)[1])

## phylo vcov by branch
for (i in seq_along(submat_size_start)) {
  condvar_branch_array[,,i] <- 
    as.matrix(
    ranef2[, -c(1, 2, 3)][
    submat_size_start[i]:submat_size_end[i]
  , submat_size_start[i]:submat_size_end[i]]
    )
}

return(condvar_branch_array)
   
}

}

#######
### Predict new respones (low level functions)
#######
lme4_titer_pred    <- function (lme4_fits, random_effects, cond_mode_sd, day, LD, body_size, brl
                                  , last_brl, species, uncertainty_list, spec_type) {
 
  fixed_effects     <- lme4_fits[["red_model_coefs"]][["Fixed_effects"]]
  rand_eff_sd       <- lme4_fits[["red_model_coefs"]][["phylo_vcov"]]
  cit_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["Citation_vcov"]]
  inf_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["unique_line_vcov"]]  
  fixed_eff_vcov    <- lme4_fits[["red_model_coefs"]][["Fixed_vcov"]]
  
  nsamps    <- uncertainty_list[["samps"]]
  pred_vals <- matrix(data = 0, nrow = 1, ncol = nsamps)
  
  ## fixed effects
    if (uncertainty_list[["titer_model.fixed_uncer"]] == FALSE) {
  ## Multivariate Normal draw from 0 variance vcov matrix (written this way to parallel code below)
  fix_eff_matrix <- mvrnorm(nsamps
    , mu    = fixed_effects
    , Sigma = matrix(data = 0, ncol = ncol(fixed_eff_vcov), nrow = nrow(fixed_eff_vcov)))
    } else {
  ## Multivariate Normal draw from the vcov of the fixed effects
  fix_eff_matrix <- mvrnorm(nsamps
    , mu    = fixed_effects
    , Sigma = fixed_eff_vcov)
    }
  
  fix_eff1 <- matrix(rep(fix_eff_matrix[, 1], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE)
  fix_eff2 <- matrix(rep(fix_eff_matrix[, 2], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) *
    matrix(data = log(day), ncol = nsamps, nrow = length(day))
  fix_eff3 <- matrix(rep(fix_eff_matrix[, 3], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) * 
    matrix(data = day, ncol = nsamps, nrow = length(day))
  fix_eff4 <- matrix(rep(fix_eff_matrix[, 4], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) * 
    matrix(data = log(body_size), ncol = nsamps, nrow = length(day)) 
  fix_eff5 <- matrix(rep(fix_eff_matrix[, 5], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) *
    matrix(data = LD, ncol = nsamps, nrow = length(day)) 
  fix_eff6 <- matrix(rep(fix_eff_matrix[, 6], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) *
    matrix(data = log(day) * log(body_size), ncol = nsamps, nrow = length(day))
  fix_eff7 <- matrix(rep(fix_eff_matrix[, 7], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) * 
    matrix(data = day * log(body_size), ncol = nsamps, nrow = length(day))
  
  ## phylogenetic random effects
    if (uncertainty_list[["titer_model.phylo_rand"]] == FALSE) {
    ## Intercept
  ran_ef1 <- replicate(nsamps, sum(random_effects[, grep("Intercept", colnames(random_effects))[1]] * brl))
    ## log(Day)
  ran_ef2 <- replicate(nsamps, sum(random_effects[, grep("log", colnames(random_effects))[1]] * brl))
    ## Day
  ran_ef3 <- replicate(nsamps, sum(random_effects[, grep("Day", colnames(random_effects))[2]] * brl))

    ## Uncertainty due to branch variation
    } else {

    ## Build an array, where the third dimension is the number of branches
   ran_ef_array <- array(data = 0, dim = c(nsamps, ncol(random_effects), length(brl)))
   
   ## random multivariate normal draws (independent for each branch)
  for (i in 1:dim(ran_ef_array)[3]) {
    ran_ef_array[,,i] <-  mvrnorm(nsamps
    , mu    = random_effects[i, ]    * brl[i]
    , Sigma = cond_mode_sd[,,i]      * brl[i]^2
      )
  }
  
   ## Sum across branches to get the [nsamps] samples for each random effect
  ran_ef_matrix <- rowSums(ran_ef_array, dim = 2) ## sum over the third dimension
  ran_ef1 <- ran_ef_matrix[, 1]
  ran_ef2 <- ran_ef_matrix[, 2]
  ran_ef3 <- ran_ef_matrix[, 3]
      
    }
  
   ran_ef1 <- matrix(data = rep(ran_ef1, length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE)
   ran_ef2 <- matrix(data = rep(ran_ef2, length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE)
   ran_ef3 <- matrix(data = rep(ran_ef3, length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE)
  
    ## Uncertainty due to tip variation. Only add this type of variation for species missing from the phylogeny
    if (uncertainty_list[["titer_model.phylo_tip"]] == TRUE & spec_type == "unknown_spec") {
      
  ## Tip variation is given by the random effect estimates for species * the last branchlength (sd given in per unit branchlength) 
    ## expectation of 0 because of the assumption of random walk
  ran_ef_matrix <- mvrnorm(nsamps
      , mu    = rep(0, nrow(rand_eff_sd))
      , Sigma = rand_eff_sd * last_brl)
  ran_ef1 <- ran_ef1 + matrix(rep(t(ran_ef_matrix[, 1]), length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)
  ran_ef2 <- ran_ef2 + matrix(rep(t(ran_ef_matrix[, 2]), length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)
  ran_ef3 <- ran_ef3 + matrix(rep(t(ran_ef_matrix[, 3]), length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)
  }
  
  ran_ef2 <- ran_ef2 * matrix(data = log(day), ncol = nsamps, nrow = length(day))
  ran_ef3 <- ran_ef3 * matrix(data = day, ncol = nsamps, nrow = length(day))
  
  ## Estiamte
    if (uncertainty_list[["titer_model.other_rand"]] == FALSE) {
  pred_vals <- exp(fix_eff1 + fix_eff2 + fix_eff3 + fix_eff4 + fix_eff5 + fix_eff6 + fix_eff7 + ran_ef1  + ran_ef2 + ran_ef3)
   } else {
  ## + Other random effects (citation and infection experiment)
  ran_ef4 <- rnorm(nsamps, 0, cit_rand_eff_sd[["sdcor"]])   
  ran_ef4 <- matrix(data = rep(ran_ef4, length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)
  ran_ef5 <- rnorm(nsamps, 0, inf_rand_eff_sd[["sdcor"]]) 
  ran_ef5 <- matrix(data = rep(ran_ef5, length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)

  pred_vals <- exp(
    fix_eff1 + fix_eff2 + fix_eff3 + fix_eff4 + fix_eff5 + fix_eff6 + fix_eff7 +
      ran_ef1  + ran_ef2 + ran_ef3  + ran_ef4 + ran_ef5)  
  
   }
  
  info_mat <- data.frame(
    species   = rep(species, length(day))
  , day       = day
  , log_dose  = LD)
  
  info_mat <- transform(info_mat, species = as.character(species))
  
  cbind(info_mat, pred_vals)
  
}
lme4_survival_pred <- function (lme4_fits, random_effects, cond_mode_sd, day, LD, body_size, brl
                                  , last_brl, species, uncertainty_list, spec_type, other_responses) {
  
## Recover the predicted titer for use in predicting survival
 ## for consistency use predicted titer for both 
names(other_responses)[1] <- "Scientific_Name"
spec_titer <- other_responses %>% 
  filter(Scientific_Name == species) %>%
  dplyr::select(-Scientific_Name, -day, -log_dose, -model)

  fixed_effects     <- lme4_fits[["red_model_coefs"]][["Fixed_effects"]]
  rand_eff_sd       <- lme4_fits[["red_model_coefs"]][["phylo_vcov"]]
  cit_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["Citation_vcov"]]
  inf_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["unique_line_vcov"]]  
  fixed_eff_vcov    <- lme4_fits[["red_model_coefs"]][["Fixed_vcov"]]
      
  nsamps    <- uncertainty_list[["samps"]]
  pred_vals <- matrix(data = 0, nrow = length(day), ncol = nsamps)
  
  ## fixed effects
    if (uncertainty_list[["survival_model.fixed_uncer"]] == FALSE) {
  ## Multivariate Normal draw from 0 variance vcov matrix (written this way to parallel code below)
  fix_eff_matrix <- mvrnorm(nsamps
    , mu    = fixed_effects
    , Sigma = matrix(data = 0, ncol = ncol(fixed_eff_vcov), nrow = nrow(fixed_eff_vcov)))
    } else {
  ## Multivariate Normal draw from the vcov of the fixed effects
  fix_eff_matrix <- mvrnorm(nsamps
    , mu    = fixed_effects
    , Sigma = fixed_eff_vcov)
    }
  
  ## multiply through fixed eff here
  fix_eff1 <- matrix(rep(fix_eff_matrix[, 1], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE)

  ## use predicted titer from the previous model
  fix_eff2 <- matrix(rep(fix_eff_matrix[, 2], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) *
     spec_titer 
  
  fix_eff3 <- matrix(rep(fix_eff_matrix[, 3], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) *
    matrix(data = log(body_size), ncol = nsamps, nrow = length(day))
  
  fix_eff4 <- matrix(rep(fix_eff_matrix[, 4], length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE) *
    matrix(data = day, ncol = nsamps, nrow = length(day))
  
    ## phylogenetic random effects
    if (uncertainty_list[["survival_model.phylo_rand"]] == FALSE) {
   
    ## Intercept   
  ran_ef1 <- replicate(nsamps, sum(random_effects[1, ] * brl))
  
    } else {
      
   ## Build an array, where the third dimension is the number of branches
   ran_ef_array <- array(data = 0, dim = c(nsamps, 1, length(brl)))
   
  ## random multivariate normal draws (independent for each branch)
  for (i in 1:dim(ran_ef_array)[3]) {
    ran_ef_array[,,i] <-  rnorm(nsamps
    , mean = random_effects[1, i]       * brl[i]
    , sd   = sqrt(cond_mode_sd[,,i])    * brl[i]^2
      )
  }
  
   ## Sum across branches to get the [nsamps] samples for each random effect
  ran_ef_matrix <- rowSums(ran_ef_array, dim = 2) ## sum over the third dimension
  ran_ef1 <- ran_ef_matrix[, 1]
  
  }
  
   ran_ef1 <- matrix(data = rep(ran_ef1, length(day)), nrow = length(day), ncol = nsamps, byrow = TRUE)
  
    ## Uncertainty due to tip variation
    if (uncertainty_list[["survival_model.phylo_tip"]] == TRUE & spec_type == "unknown_spec") {
    ## This whole function is already very model specific, so while numbers here are not optimal they will do for now
      ## sqrt(rand_eff_sd) --> extracted as vcov, but only 1 random effect so just var returned
  ran_ef_matrix <- rnorm(nsamps, 0, sqrt(rand_eff_sd) * last_brl)
  ran_ef1 <- ran_ef1 + matrix(rep(ran_ef_matrix, length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)
    }

    ## Estiamte
    if (uncertainty_list[["survival_model.other_rand"]] == FALSE) {
      
    ## written to be explicit
  pred_vals <- 1 - (1 - exp(-exp(fix_eff1 + fix_eff2 + fix_eff3 + fix_eff4 + ran_ef1)))
  
    } else {
  ## + Other random effects (citation and infection experiment)
  ran_ef2 <- rnorm(nsamps, 0, cit_rand_eff_sd[["sdcor"]])   
  ran_ef2 <- matrix(data = rep(ran_ef2, length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)
  ran_ef3 <- rnorm(nsamps, 0, inf_rand_eff_sd[["sdcor"]]) 
  ran_ef3 <- matrix(data = rep(ran_ef3, length(day)), ncol = nsamps, nrow = length(day), byrow = TRUE)  
     
    ## written to be explicit 
  pred_vals     <- 1 - (1 - exp(-exp(fix_eff1 + fix_eff2 + fix_eff3 + fix_eff4 + ran_ef1 + ran_ef2 + ran_ef3))) 
    }
  
  ## Calculate cumulative survival
  pred_vals     <- apply(pred_vals, 2, cumprod)
  
  info_mat <- data.frame(
    species   = rep(species, length(day))
  , day       = day
  , log_dose  = LD
  , titer     = rowMeans(spec_titer)) ## Just a placeholder because all 1000 samples are used
  info_mat <- transform(info_mat, species = as.character(species))
  
  cbind(info_mat, pred_vals)
}
lme4_biting_pred   <- function (lme4_fits, random_effects, cond_mode_sd, body_size, brl, last_brl
                                  , species, uncertainty_list, spec_type) {
  
  ## Keeping code consistent with the other functions of (nrow = length(day) with "1")
  nsamps    <- uncertainty_list[["samps"]]
  pred_vals <- matrix(data = 0, nrow = 1, ncol = nsamps)
  
  fixed_effects     <- lme4_fits[["red_model_coefs"]][["Fixed_effects"]]
  rand_eff_sd       <- lme4_fits[["red_model_coefs"]][["phylo_vcov"]]
  cit_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["Citation_vcov"]]
  inf_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["unique_line_vcov"]]  
  fixed_eff_vcov    <- lme4_fits[["red_model_coefs"]][["Fixed_vcov"]] 
  
  ## fixed effects
    if (uncertainty_list[["bite_model.fixed_uncer"]] == FALSE) {
  ## Normal draw from 0 variance vcov matrix (written this way to parallel code below)
  fix_eff_matrix <- rnorm(nsamps, fixed_effects, 0)
    } else {
  ## Normal draw from the vcov of the fixed effects
  fix_eff_matrix <- rnorm(nsamps, fixed_effects, fixed_eff_vcov[1, 1])
    }
  
  ## multiply through fixed eff here ***
  fix_eff1 <- matrix(rep(fix_eff_matrix, 1), nrow = 1, ncol = nsamps, byrow = TRUE)
  
    ## phylogenetic random effects
    if (uncertainty_list[["bite_model.phylo_rand"]] == FALSE) {
  ran_ef1 <- replicate(nsamps, sum(random_effects[1, ] * brl))
    } else {
      
  ## Build an array, where the third dimension is the number of branches
   ran_ef_array <- array(data = 0, dim = c(nsamps, 1, length(brl)))
   
  ## random multivariate normal draws (independent for each branch)
  for (i in 1:dim(ran_ef_array)[3]) {
    ran_ef_array[,,i] <-  rnorm(nsamps
    , mean = random_effects[1, i]       * brl[i]
    , sd   = sqrt(cond_mode_sd[,,i])    * brl[i]^2
      )
  }
  
  ## Sum across branches to get the [nsamps] samples for each random effect
  ran_ef_matrix <- rowSums(ran_ef_array, dim = 2) ## sum over the third dimension
  ran_ef1 <- ran_ef_matrix[, 1]
  
  }
  
  ran_ef1 <- matrix(data = rep(ran_ef1, 1), nrow = 1, ncol = nsamps, byrow = TRUE)
  
  ## Uncertainty due to tip variation
    if (uncertainty_list[["bite_model.phylo_tip"]] == TRUE & spec_type == "unknown_spec") {
  ## This whole function is already very model specific, so while numbers here are not optimal they will do for now
  ran_ef_matrix <- rnorm(nsamps, 0, sqrt(rand_eff_sd) * last_brl)
  ran_ef1 <- ran_ef1 + matrix(rep(ran_ef_matrix, 1), ncol = nsamps, nrow = 1, byrow = TRUE)
    }
  
  ## Estiamte
  pred_vals <- exp(fix_eff1 + ran_ef1)

  info_mat <- data.frame(
    species   = species)
  info_mat <- transform(info_mat, species = as.character(species))
  
  cbind(info_mat, pred_vals)
}
lme4_detect_pred   <- function (lme4_fits, random_effects, cond_mode_sd, body_size, brl, last_brl
                                  , species, uncertainty_list, spec_type) {
 
  ## Keeping code consistent with the other functions of (nrow = length(day) with "1")
  nsamps <- uncertainty_list[["samps"]]
  pred_vals <- matrix(data = 0, nrow = 1, ncol = nsamps)
  
  fixed_effects     <- lme4_fits[["red_model_coefs"]][["Fixed_effects"]]
  rand_eff_sd       <- lme4_fits[["red_model_coefs"]][["phylo_vcov"]]
  cit_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["Citation_vcov"]]
  inf_rand_eff_sd   <- lme4_fits[["red_model_coefs"]][["unique_line_vcov"]]  
  fixed_eff_vcov    <- lme4_fits[["red_model_coefs"]][["Fixed_vcov"]]
  
  ## fixed effects
      if (uncertainty_list[["detect_model.fixed_uncer"]] == FALSE) {
  ## Multivariate Normal draw from 0 variance vcov matrix (written this way to parallel code below)
  fix_eff_matrix <- mvrnorm(nsamps
    , mu    = fixed_effects
    , Sigma = matrix(data = 0, ncol = ncol(fixed_eff_vcov), nrow = nrow(fixed_eff_vcov)))
    } else {
  ## Multivariate Normal draw from the vcov of the fixed effects
  fix_eff_matrix <- mvrnorm(nsamps
    , mu    = fixed_effects
    , Sigma = fixed_eff_vcov)
    }
  
  fix_eff1 <- matrix(rep(fix_eff_matrix[, 1], 1), nrow = 1, ncol = nsamps, byrow = TRUE)
  fix_eff2 <- matrix(rep(fix_eff_matrix[, 2], 1), nrow = 1, ncol = nsamps, byrow = TRUE) *
    matrix(data = log(body_size), ncol = nsamps, nrow = 1)
  
  ## phylogenetic random effects
      if (uncertainty_list[["bite_model.phylo_rand"]] == FALSE) {
  ran_ef1 <- replicate(nsamps, sum(random_effects[1, ] * brl))
    } else {
      
  ## Build an array, where the third dimension is the number of branches
   ran_ef_array <- array(data = 0, dim = c(nsamps, 1, length(brl)))
   
  ## random multivariate normal draws (independent for each branch)
  for (i in 1:dim(ran_ef_array)[3]) {
    ran_ef_array[,,i] <-  rnorm(nsamps
    , mean = random_effects[1, i]       * brl[i]
    , sd   = sqrt(cond_mode_sd[,,i])    * brl[i]^2
      )
  }
  
  ## Sum across branches to get the [nsamps] samples for each random effect
  ran_ef_matrix <- rowSums(ran_ef_array, dim = 2) ## sum over the third dimension
  ran_ef1 <- ran_ef_matrix[, 1]
  
    }
  
  ran_ef1 <- matrix(data = rep(ran_ef1, 1), nrow = 1, ncol = nsamps, byrow = TRUE)
  
  ## Uncertainty due to tip variation
    if (uncertainty_list[["bite_model.phylo_tip"]] == TRUE & spec_type == "unknown_spec") {
  ## This whole function is already very model specific, so while numbers here are not optimal they will do for now
  ran_ef_matrix <- rnorm(nsamps, 0, sqrt(rand_eff_sd) * last_brl)
  ran_ef1 <- ran_ef1 + matrix(rep(ran_ef_matrix, 1), ncol = nsamps, nrow = 1, byrow = TRUE)
    }
  
    ## Estiamte
  pred_vals <- exp(fix_eff1 + fix_eff2 + ran_ef1)
  
  info_mat <- data.frame(
    species   = species
  , body_size = body_size)
  info_mat <- transform(info_mat, species = as.character(species))
  
  cbind(info_mat, pred_vals)
}

#######
### Wrappers at various levels (mid level functions)
#######
lme4_estimator_logcial         <- function (lme4_fits, random_effects, cond_mode_sd, day, LD, body_size, brl, last_brl
                                           , species, which_model, uncertainty_list, spec_type
                                           , other_responses) {
    ## call estimator function
if (which_model == "titer") {

  lme4_titer_pred(
    lme4_fits         = lme4_fits
  , random_effects    = random_effects
  , cond_mode_sd      = cond_mode_sd
  , day               = day
  , LD                = LD
  , body_size         = body_size
  , brl               = brl
  , last_brl          = last_brl  
  , species           = species
  , uncertainty_list  = uncertainty_list
  , spec_type         = spec_type)
  
} else if (which_model == "survival") {
  
  lme4_survival_pred(
    lme4_fits         = lme4_fits
  , random_effects    = random_effects
  , cond_mode_sd      = cond_mode_sd
  , day               = day
  , LD                = LD
  , body_size         = body_size
  , brl               = brl
  , last_brl          = last_brl  
  , species           = species
  , uncertainty_list  = uncertainty_list
  , spec_type         = spec_type
  , other_responses   = other_responses)
  
} else if (which_model == "bite") {
  
  lme4_biting_pred(
    lme4_fits         = lme4_fits
  , random_effects    = random_effects
  , cond_mode_sd      = cond_mode_sd
  , body_size         = body_size
  , brl               = brl
  , last_brl          = last_brl  
  , species           = species
  , uncertainty_list  = uncertainty_list
  , spec_type         = spec_type)
  
} else if (which_model == "detect") {
  
  lme4_detect_pred(
    lme4_fits         = lme4_fits
  , random_effects    = random_effects
  , cond_mode_sd      = cond_mode_sd
  , body_size         = body_size
  , brl               = brl
  , last_brl          = last_brl  
  , species           = species
  , uncertainty_list  = uncertainty_list
  , spec_type         = spec_type)
  
}
  
}
lme4_model_logcial          <- function (data, phylo, phyloZ, nsp, which_model, uncertainty_list) {
 
  if (which_model == "titer") {
    
lme4_titer_model_fit(
     data     = data
  ,  phylo    = phylo
  ,  phyloZ   = phyloZ
  ,  nsp      = nsp)
    
} else if (which_model == "survival") {
  
lme4_survival_model_fit(
     data     = data
  ,  phylo    = phylo
  ,  phyloZ   = phyloZ
  ,  nsp      = nsp) 
  
} else if (which_model == "bite") {
  
lme4_biting_model_fit(
     data             = data
  ,  phylo            = phylo
  ,  phyloZ           = phyloZ
  ,  nsp              = nsp
  ,  uncertainty_list = uncertainty_list)

} else if (which_model == "detect") {
  
lme4_detect_model_fit(
     data     = data
  ,  phylo    = phylo
  ,  phyloZ   = phyloZ
  ,  nsp      = nsp)
  
}
  
}
lme4_rand_ef_vec_logical    <- function (model_coefs, rand_eff_est, rand_ef_id, rand_ef_lengths
                                           , names_vec, which_model, lme4_fit, phylo_data, is_sd) {
  
if (which_model == "titer") {
  
lme4_titer_rand_ef_vec(
    model_coefs     = model_coefs
  , rand_eff_est    = rand_eff_est
  , rand_ef_id      = rand_ef_id
  , rand_ef_lengths = rand_ef_lengths
  , names_vec       = names_vec
  , lme4_fit        = lme4_fit
  , phylo_data      = phylo_data
  , is_sd           = is_sd)
  
} else if (which_model == "survival") {
  
lme4_survival_rand_ef_vec(
    model_coefs     = model_coefs
  , rand_eff_est    = rand_eff_est
  , rand_ef_id      = rand_ef_id
  , rand_ef_lengths = rand_ef_lengths
  , names_vec       = names_vec
  , lme4_fit        = lme4_fit
  , phylo_data      = phylo_data
  , is_sd           = is_sd)
  
} else if (which_model == "bite") {
  
lme4_biting_rand_ef_vec(
    model_coefs     = model_coefs
  , rand_eff_est    = rand_eff_est
  , rand_ef_id      = rand_ef_id
  , rand_ef_lengths = rand_ef_lengths
  , names_vec       = names_vec
  , lme4_fit        = lme4_fit
  , phylo_data      = phylo_data
  , is_sd           = is_sd)
  
} else if (which_model == "detect") {
  
lme4_detect_rand_ef_vec(
    model_coefs     = model_coefs
  , rand_eff_est    = rand_eff_est
  , rand_ef_id      = rand_ef_id
  , rand_ef_lengths = rand_ef_lengths
  , names_vec       = names_vec
  , lme4_fit        = lme4_fit
  , phylo_data      = phylo_data
  , is_sd           = is_sd)
  
}
  
}
## Logical to load correct previously saved estimated responses
response_load_logcial       <- function (which_model, uncertainty_list) {
  
    ## Build a name to save the results as
  uncer_name      <- vector("character", length(uncertainty_list))
  for (i in 1:length(uncertainty_list)) {
    uncer_name[i] <- paste(names(uncertainty_list)[i], as.character(uncertainty_list[[i]]), sep = "_")
  }

  if (which_model == "titer") {
    
    ## remove the uncertainty types I don't want to save
    
    ## If uncertainty_list[[XXXX_model]] == TRUE, check for the full name to make sure the correct combination
      ## of uncertainty is found. 
    ## If uncertainty_list[[XXXX_model]] == FALSE, then there is no uncertainty in that model, and the rest of 
      ## the options don't matter
  
    if (uncertainty_list[["titer_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["titer_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["titer_model.other_rand"]]  == TRUE
      | uncertainty_list[["titer_model.phylo_tip"]]   == TRUE) {
  
    needed_uncer     <- grep(paste(c("titer_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("titer_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")
    
    } else {
      
    needed_uncer     <- grep("titer_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("titer_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")      
      
    }
    
} else if (which_model == "survival") {
  
    if (uncertainty_list[["survival_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["survival_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["survival_model.other_rand"]]  == TRUE
      | uncertainty_list[["survival_model.phylo_tip"]]   == TRUE) {
  
    needed_uncer     <- grep(paste(c("survival_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("survival_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")
    
    } else {
      
    needed_uncer     <- grep("survival_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("survival_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")      
      
    }
    
} else if (which_model == "bite") {
  
    if (uncertainty_list[["bite_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["bite_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["bite_model.phylo_tip"]]   == TRUE) {
  
    needed_uncer     <- grep(paste(c("bite_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("bite_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")
    
    } else {
      
    needed_uncer     <- grep("bite_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("bite_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")      
      
    }
  
} else if (which_model == "detect") {
  
    if (uncertainty_list[["detect_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["detect_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["detect_model.phylo_tip"]]   == TRUE)  {
  
    needed_uncer     <- grep(paste(c("detect_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("detect_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")
    
    } else {
      
    needed_uncer     <- grep("detect_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("detect_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")
      
    }
  
}
  
  if (file.exists(full.file.name))  {
    prev_fits <- read.csv(full.file.name, sep = "")
    prev_fits <- transform(prev_fits, species = as.character(species))
  }
    
}
## Logical to check for previously saved estimated responses in loaded data frame from above function
response_check_logcial      <- function (prev_responses, species) {
  
  resp_return <- prev_responses[prev_responses[["species"]] == species, ]
  transform(resp_return, species = as.character(species))
    
}
## Logical to write over correct previously saved estimated responses
response_write_logcial      <- function (which_model, new_species_estimate, prev_responses, uncertainty_list) {
  
  ## Build a name to save the results as
  uncer_name      <- vector("character", length(uncertainty_list))
  for (i in 1:length(uncertainty_list)) {
    uncer_name[i] <- paste(names(uncertainty_list)[i], as.character(uncertainty_list[[i]]), sep = "_")
  }
  
  if (which_model == "titer") {
    
    if (uncertainty_list[["titer_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["titer_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["titer_model.other_rand"]]  == TRUE
      | uncertainty_list[["titer_model.phylo_tip"]]   == TRUE) {
    
    ## Select the uncertainty types that I want to save
    needed_uncer     <- grep(paste(c("titer_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("titer_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")
    
    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1, 2)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }
    
   } else {
   
    needed_uncer     <- grep("titer_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("titer_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1, 2)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }
     
   }

} else if (which_model == "survival") {
  
    if (uncertainty_list[["survival_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["survival_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["survival_model.other_rand"]]  == TRUE
      | uncertainty_list[["survival_model.phylo_tip"]]   == TRUE) {
  
    needed_uncer     <- grep(paste(c("survival_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("survival_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1, 2)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }
    
    } else {
    
    needed_uncer     <- grep("survival_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("survival_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1, 2)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }
      
    }

} else if (which_model == "bite") {
  
    if (uncertainty_list[["bite_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["bite_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["bite_model.phylo_tip"]]   == TRUE) {
  
    needed_uncer     <- grep(paste(c("bite_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("bite_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }
    
    } else {
      
    needed_uncer     <- grep("bite_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("bite_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }
       
    }

} else if (which_model == "detect") {
  
    if (uncertainty_list[["detect_model.fixed_uncer"]] == TRUE 
      | uncertainty_list[["detect_model.phylo_rand"]]  == TRUE
      | uncertainty_list[["detect_model.phylo_tip"]]   == TRUE) {
  
    needed_uncer     <- grep(paste(c("detect_model"), collapse = "|"), uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("detect_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }

    } else {
  
    needed_uncer     <- grep("detect_model", uncer_name)
    uncer_name       <- as.character(glue::collapse(uncer_name[c(needed_uncer)], sep = "_"))
    file.name        <- paste("detect_resp", uncer_name, sep = "_")
    full.file.name   <- paste("saved_fits/", file.name, sep = "")

    ## Add new fits to the old fits
    if (!is.null(prev_responses)) {
    all_fits <- rbind(new_species_estimate, prev_responses)
    all_fits <- all_fits[!duplicated(all_fits[, c(1)]), ]
    write.table(all_fits, file = full.file.name)
    } else {
    write.table(new_species_estimate, file = full.file.name)  
    }  
      
    }
    
}
    
}
## Subset to only the species that I have data for to run leave-one-out validation
leave_one_out_setup_logical <- function (data, phylo) {

## species that I have data for
leave_one_spec     <- unique(data[["Scientific_Name"]])
leave_one_spec_vec <- numeric(length(leave_one_spec))

## tip # associated with each species
for (z in seq_along(leave_one_spec)) {
 leave_one_spec_vec[z] <- which(sing_tree[["tip.label"]] == as.character(leave_one_spec[z]))
}

## The tips I have and the tips I need to estiamte for are the same for leave-one-out validation
specs_to_est <- list(
  needed_tips  = leave_one_spec_vec
, have_tips    = leave_one_spec_vec
, missing_spec = seq(1, length(leave_one_spec_vec), by = 1)
)

return(specs_to_est)

}

#######
### Extract parameters and make predictions (mid level functions)
#######
lme4_model_run_extract         <- function (data.full, phylo_data, which_model, fit_full, fit_red
                                              , lme4_fit_full, spec_name, outside, uncertainty_list) {
  
## fit whichever models (full and reduced) the specific call asks for
if (fit_full == TRUE) {
  
lme4_fit_full <- lme4_model_logcial(
  data             = data.full
, phylo            = phylo_data[["full_phylo"]]
, phyloZ           = phylo_data[["full_Z_mat"]]
, nsp              = nrow(phylo_data[["full_Z_mat"]])
, which_model      = which_model
, uncertainty_list = uncertainty_list)

full_model_coefs <- lme4_extract_est(
  data             = data.full
, lme4_fit         = lme4_fit_full
, phylo            = phylo_data[["full_phylo"]]
, phyloZ           = phylo_data[["full_Z_mat"]]
, which_model      = which_model)
}

if (fit_red == TRUE) {
  
  ## If leave-one-out is being run, subset the data
if (outside == FALSE) {
  data.use <- droplevels(subset(data.full, Scientific_Name != spec_name))
} else {
  data.use <- data.full
}

lme4_fit_red <- lme4_model_logcial(
  data             = data.use
, phylo            = phylo_data[["reduced_phylo"]]
, phyloZ           = phylo_data[["reduced_Z_mat"]]
, nsp              = nrow(phylo_data[["reduced_Z_mat"]])
, which_model      = which_model
, uncertainty_list = uncertainty_list)

red_model_coefs <- lme4_extract_est(
     data        = data.use   
  ,  lme4_fit    = lme4_fit_red
  ,  phylo       = phylo_data[["reduced_phylo"]]
  ,  phyloZ      = phylo_data[["reduced_Z_mat"]]
  ,  which_model = which_model)

}

## empty storage if model wasn't run
if (!exists("lme4fit_full")) lme4fit_full         <- NULL
if (!exists("lme4fit_red")) lme4fit_red           <- NULL
if (!exists("full_model_coefs")) full_model_coefs <- NULL
if (!exists("red_model_coefs")) red_model_coefs   <- NULL

return(list(
  lme4_fit_full    = lme4_fit_full
, lme4_fit_red     = lme4_fit_red
, full_model_coefs = full_model_coefs
, red_model_coefs  = red_model_coefs
))

}
lme4_extract_est               <- function (data, lme4_fit, phylo, phyloZ, which_model) {
  
## random effect estimates
rand_eff_est    <- getME(lme4_fit, c("b"))@x

## break up the random effects vector into the estimates for each random effect by looking at the order of the random effects
rand_ef_lengths <- numeric(length = length(unlist(lme4_fit@cnms)))

## expand lme4fit@cnms to repeat sci-name (special case of correlated random effects for titer model)
#ran_ef_terms <- apply(matrix(unique(rand_ef_lengths)), 1, function(x) length(which(rand_ef_lengths == x)))
ran_ef_terms <- unlist(lapply(lme4_fit@cnms, function (x) length(x)))
names_vec    <- rep(names(lme4_fit@cnms), ran_ef_terms)

## Store number of phylogenetic random effects (estimates + sd)
total_ran_efs    <- length(names_vec)

## determine the lengths of each random effect
for (i in seq_along(names_vec)) {
  rand_ef_lengths[i] <- ifelse(
    names_vec[i] == "Scientific_Name"
  , ncol(phyloZ)  ## number of branches in the phylogeny
  , length(unique(data[[names_vec[i]]])) ## length of the unique values in the data
  )
}

## Separate out the random effects
   ## +4 for fixed effects, summary matrix of model coefs, random ef sds and fixed ef vocv
model_coefs <- vector("list", length = total_ran_efs + 4 + length(ran_ef_terms) + 1)
rand_ef_id  <- rbind(names_vec, unlist(lme4_fit@cnms))

## organize the random effect conditional modes appropriately. Call a logical sorting function
  ## depending on the model being run
model_coefs <- lme4_rand_ef_vec(
    model_coefs     = model_coefs
  , rand_eff_est    = rand_eff_est
  , rand_ef_id      = rand_ef_id
  , rand_ef_lengths = rand_ef_lengths
  , names_vec       = names_vec
  , lme4_fit        = lme4_fit
  , phylo_data      = phyloZ
  , is_sd           = FALSE) 

## Extract sd on the conditional modes of the random effects
  ## recover only those for the phylogeny, the rest are not being used
## https://github.com/lme4/lme4/issues/148
#cond_cov_mat <- ranef(lme4_fit, condVar = TRUE)
#cond_cov_mat <- attr(cond_cov_mat[["phylo"]],"postVar")
cond_cov_mat   <- lme4:::condVar(lme4_fit)
phylo_coefs_sd <- vector("list", length = total_ran_efs)
  
phylo_coefs_sd <- lme4_rand_ef_vec(
    model_coefs     = phylo_coefs_sd
  , rand_eff_est    = cond_cov_mat
  , rand_ef_id      = rand_ef_id
  , rand_ef_lengths = rand_ef_lengths
  , names_vec       = names_vec
  , lme4_fit        = lme4_fit
  , phylo_data      = phyloZ
  , is_sd           = TRUE) 

## For the second to last entry combine the phylo random effect into a matrix
phylo_ran_ef <- grep("Scientific_Name", names(model_coefs))

model_coefs[[total_ran_efs + length(ran_ef_terms) + 2]] <- matrix(
  ncol     = length(phylo_ran_ef)
, nrow     = ncol(phyloZ)
, data     = unlist(model_coefs[phylo_ran_ef])
, dimnames = list(NULL, names(model_coefs[phylo_ran_ef])))

## For the last entry include the vector of fixed effects
model_coefs[[total_ran_efs + length(ran_ef_terms) + 3]] <- getME(lme4_fit, c("beta"))

## fixed effects vcov 
model_coefs[[total_ran_efs + length(ran_ef_terms) + 4]] <- summary(lme4_fit)[["vcov"]]

## vcov of the conditional random effects
model_coefs[[total_ran_efs + length(ran_ef_terms) + 5]] <- phylo_coefs_sd

## Fill in names 
names(model_coefs)[
  (total_ran_efs + length(ran_ef_terms) + 2):
  (total_ran_efs + length(ran_ef_terms) + 5)] <- c(
    "Phylogenetic_random_effects"
  , "Fixed_effects"
  , "Fixed_vcov"
  , "Conditional_covariance"
    )

return(model_coefs)

}
  ## species without data
estiamte_response_unknown_spec <- function (phylo_data, lme4_fits, data.full
                                              , species, species_body_size
                                              , which_model, uncertainty_list
                                              , other_responses) {
## set up branches and coefficients needed
  ## Predict for reduced model by selecting random effect for species | spec Z vec
model_coefs_needed <- lme4_fits[["red_model_coefs"]][["Phylogenetic_random_effects"]][which(phylo_data[["species_vector"]] > 0), ]

## Conditional vcov
cond_mode_sd  <- lme4_fits[["red_model_coefs"]][["Conditional_covariance"]][,,which(phylo_data[["species_vector"]] > 0)]

ran_ef_terms  <- unlist(lapply(lme4_fits$lme4_fit_red@cnms, function (x) length(x)))
num_phylo_ran <- ran_ef_terms[grep("Scientific_Name", names(ran_ef_terms))]

## Needed for species that are an out group
## check if not in vcov structure (no correlated random effects or an outgroup (one branch))
if (length(dim(cond_mode_sd)) != 3) {
  if (length(dim(cond_mode_sd)) != 2) {
   cond_mode_sd <- array(data = cond_mode_sd, dim = c(1, 1, length(cond_mode_sd)))  
  } else {
   cond_mode_sd <- array(data = cond_mode_sd, dim = c(num_phylo_ran, num_phylo_ran, 1))
  }
}

## These are the branchlengths that need to be multiplied by the estimates. Note, for a normal random effect these are 
 ## indicators, but here they are != 1 repeat them equal to the number of random effects
brl_for_spec <- phylo_data[["species_vector"]][which(phylo_data[["species_vector"]] != 0)]
  
## if a vector (single branch) create a matrix for the calculation
if (class(model_coefs_needed) == "numeric") {
  temp_colnames                <- names(model_coefs_needed)
  model_coefs_needed           <- matrix(model_coefs_needed, ncol = length(model_coefs_needed)
    , nrow = 1, byrow = TRUE)
  colnames(model_coefs_needed) <- temp_colnames 
}
  
## Set up data for predictions (wont need these for all response variables so have a logical here)
if (which_model == "titer" | which_model == "survival") {
pred_day <- seq(1, 8, by = 1)
#pred_LD  <- rep(mean(data.full[["Log_Dose"]]), length(pred_day))
## Mean dose given by C. tarsalis
pred_LD  <- rep(6.1, length(pred_day))
} else {
pred_day <- NULL
pred_LD  <- NULL
}
body_size   <- community_spec_bs[community_spec_bs[["species"]] == as.character(species), ][["body_size_s"]]
if (length(body_size) == 0) {
  body_size <- mean(community_spec_bs[["body_size_s"]])
}

pred_vals <- lme4_estimator_logcial(
  lme4_fits         = lme4_fits
, random_effects    = model_coefs_needed
, cond_mode_sd      = cond_mode_sd
, day               = pred_day
, LD                = pred_LD
, body_size         = body_size
, brl               = brl_for_spec
, last_brl          = phylo_data[["last_brl"]]
, species           = species
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, spec_type         = "unknown_spec"
, other_responses   = other_responses)

## stupid factors
pred_vals <- transform(pred_vals, species = as.character(species))

## return a data frame with the first 4 columns as the "info" columns
pred_vals_r <- cbind(
  pred_vals[, 1:(ncol(pred_vals) - uncertainty_list[["samps"]])]
, data.frame(model = rep("unknown_response", nrow(pred_vals)))
, pred_vals[, (1 + ncol(pred_vals) - uncertainty_list[["samps"]]):
    ((ncol(pred_vals) - uncertainty_list[["samps"]]) + uncertainty_list[["samps"]])])
names(pred_vals_r)[1] <- "species"

return(pred_vals_r)

}
  ## species with data
estimate_response_known_spec   <- function (phylo_data, lme4_fits, data.full 
                                              , species, species_body_size
                                              , which_model, outside, uncertainty_list
                                              , other_responses) {
  
## different data is used whether or not leave-one-out or outside est is being run
if (outside == TRUE) {
 
## (1) Predict for reduced model by selecting random effect for species | spec Z vec
  model_coefs_needed <- lme4_fits[["red_model_coefs"]][["Phylogenetic_random_effects"]][
    which(phylo_data[["reduced_Z_mat"]][which(phylo_data[["reduced_Z_mat"]]@Dimnames[[1]] == species), ] != 0), ]

## Conditional vcov
cond_mode_sd   <- lme4_fits[["red_model_coefs"]][["Conditional_covariance"]][,,
  which(phylo_data[["reduced_Z_mat"]][which(phylo_data[["reduced_Z_mat"]]@Dimnames[[1]] == species), ] != 0)]

ran_ef_terms  <- unlist(lapply(lme4_fits$lme4_fit_red@cnms, function (x) length(x)))
num_phylo_ran <- ran_ef_terms[grep("Scientific_Name", names(ran_ef_terms))]

## Needed for species that are an out group
## check if not in vcov structure (no correlated random effects or an outgroup (one branch))
if (length(dim(cond_mode_sd)) != 3) {
  if (length(dim(cond_mode_sd)) != 2) {
   cond_mode_sd <- array(data = cond_mode_sd, dim = c(1, 1, length(cond_mode_sd)))  
  } else {
   cond_mode_sd <- array(data = cond_mode_sd, dim = c(num_phylo_ran, num_phylo_ran, 1))
  }
}
  
  ## These are the branchlengths that need to be multiplied by the estimates. Note, for a normal
    ## random effect these are indicators, but here they are != 1; repeat them equal to the number of random effects
  brl_for_spec <- phylo_data[["reduced_Z_mat"]][which(phylo_data[["reduced_Z_mat"]]@Dimnames[[1]] == species)
    , which(phylo_data[["reduced_Z_mat"]][which(phylo_data[["reduced_Z_mat"]]@Dimnames[[1]] == species), ] != 0)]
  
## if a vector (single branch) create a matrix for the calculation
if (class(model_coefs_needed) == "numeric") {
  temp_colnames                <- names(model_coefs_needed)
  model_coefs_needed           <- matrix(model_coefs_needed, ncol = length(model_coefs_needed)
    , nrow = 1, byrow = TRUE)
  colnames(model_coefs_needed) <- temp_colnames
}
    
  } else {
  
  mod_spec_Z         <- phylo_data[["full_Z_mat"]][which(phylo_data[["full_Z_mat"]]@Dimnames[[1]] == species), ]
  model_coefs_needed <- lme4_fits[["full_model_coefs"]][["Phylogenetic_random_effects"]][which(mod_spec_Z > 0), ]
  brl_for_spec       <- mod_spec_Z[which(mod_spec_Z != 0)]    
      
  }

## Set up data for predictions (wont need these for all response variables so have a logical here)
if (which_model == "titer" | which_model == "survival") {
pred_day <- seq(1, 8, by = 1)
#pred_LD  <- rep(mean(data.full[["Log_Dose"]]), length(pred_day))
## Mean dose given by C. tarsalis
pred_LD  <- rep(6.1, length(pred_day))
} else {
pred_day <- NULL
pred_LD  <- NULL
}
body_size   <- community_spec_bs[community_spec_bs[["species"]] == as.character(species), ][["body_size_s"]]
if (length(body_size) == 0) {
  body_size <- mean(community_spec_bs[["body_size_s"]])
}
  
  ## call estimator function
pred_vals <- lme4_estimator_logcial(
  lme4_fits         = lme4_fits
, random_effects    = model_coefs_needed
, cond_mode_sd      = cond_mode_sd
, day               = pred_day
, LD                = pred_LD
, body_size         = body_size
, brl               = brl_for_spec
, last_brl          = phylo_data[["last_brl"]]
, species           = species
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, spec_type         = "known_spec"
, other_responses   = other_responses)

## stupid factors
pred_vals <- transform(pred_vals, species = as.character(species))

## return a data frame with the first 4 columns as the "info" columns
pred_vals_r           <- cbind(
  pred_vals[, 1:(ncol(pred_vals) - uncertainty_list[["samps"]])]
, data.frame(model = rep("known_response", nrow(pred_vals)))
, pred_vals[, (1 + ncol(pred_vals) - uncertainty_list[["samps"]]):
    ((ncol(pred_vals) - uncertainty_list[["samps"]]) + uncertainty_list[["samps"]])])
names(pred_vals_r)[1] <- "species"

return(pred_vals_r)

}
  ## species absent from the phylogeny
estiamte_response_missing_spec <- function (phylo_data, lme4_fits, data.full 
                                              , species, species_body_size
                                              , which_model, uncertainty_list
                                              , other_responses) {

## Set up data for predictions (wont need these for all response variables so have a logical here)
if (which_model == "titer" | which_model == "survival") {
pred_day <- seq(1, 8, by = 1)
#pred_LD  <- rep(mean(data.full[["Log_Dose"]]), length(pred_day))
## Mean dose given by C. tarsalis
pred_LD  <- rep(6.1, length(pred_day))
pred_bs  <- rep(species_body_size, length(pred_day)) 
} else {
pred_day <- NULL
pred_LD  <- NULL
pred_bs  <- species_body_size
}

ran_ef_terms       <- unlist(lapply(lme4_fits$lme4_fit_red@cnms, function (x) length(x)))
use.random_effects <- matrix(data = 0 , nrow = 1, ncol = ncol(lme4_fits$red_model_coefs$Phylogenetic_random_effects))
num_phylo_ran      <- ran_ef_terms[grep("Scientific_Name", names(ran_ef_terms))]

## single logical needed here for missing spec because of how many random effects are in each model
cond_mode_sd       <- array(data = 0, dim = c(num_phylo_ran, num_phylo_ran, 1))

pred_vals <- lme4_estimator_logcial(
  lme4_fits         = lme4_fits  
, random_effects    = use.random_effects
, cond_mode_sd      = cond_mode_sd 
, day               = pred_day
, LD                = pred_LD
, body_size         = pred_bs
, brl               = 0
, last_brl          = 0
, species           = species
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, spec_type         = "missing_spec"
, other_responses   = other_responses)

## stupid factors
pred_vals <- transform(pred_vals, species = as.factor(species))

## return a data frame with the first 4 columns as the "info" columns
pred_vals_r           <- cbind(
  pred_vals[, 1:(ncol(pred_vals) - uncertainty_list[["samps"]])]
, data.frame(model = rep("missing_response", nrow(pred_vals)))
, pred_vals[, (1 + ncol(pred_vals) - uncertainty_list[["samps"]]):
    ((ncol(pred_vals) - uncertainty_list[["samps"]]) + uncertainty_list[["samps"]])])
names(pred_vals_r)[1] <- "species"

return(pred_vals_r)

} 

#######
### Top level wrapper (top level function)
#######
outside_spec_est               <- function (phylo, spec_with_dat, have_tips, need_tips
                                              , missing_spec, which_species_no_phylo, data.full
                                              , community_spec_bs, which_model, outside
                                              , write_resp, use_saved_responses, uncertainty_list
                                              , other_responses) {

## Set up species list and timing record
full_species <- seq(1, length(phylo[["tip.label"]]))
time_check   <- numeric(length(missing_spec))

### Check to see if a species response has previously been run and stored. 
if (use_saved_responses == TRUE) {
  ## load them
prev_responses <- response_load_logcial(
  which_model      = which_model
, uncertainty_list = uncertainty_list)

if (!is.null(prev_responses)) print("Responses Loaded")
} else {
  ## store previous responses as NULL
  prev_responses <- NULL
}

## loop over all species that I need estimates for 
for (i in 1:length(missing_spec)) {
#for (i in 1:3) {  
  
time_check[i] <- system.time({
 
## Adding a single species to the species for which I have data
temp_spec <- c(have_tips, need_tips[missing_spec[i]])
temp_phy  <- drop.tip(phylo, tip = full_species[-temp_spec])
spec_name <- phylo[["tip.label"]][need_tips[missing_spec[i]]]

  ## For run i, check within the loaded data frame for the species that is about to be calculated
  if (!is.null(prev_responses)) {
   saved_spec_res <- response_check_logcial(
     prev_responses = prev_responses
  ,  species        = spec_name
   )
   ## If a species response wasn't found, remove this saved variable
   if (nrow(saved_spec_res) == 0) {
     rm(saved_spec_res)
   }
  }

## only fit the model the first time (model with all species for which I have data
  ## only needs to be run a single time)
if (i == 1) {
  
z_mats_phylos <- recov_spec_z(
  full_phylo = temp_phy
, species    = spec_name) 

lme4_fits <- lme4_model_run_extract(
  data.full        = data.full
, phylo_data       = z_mats_phylos
, which_model      = which_model
, fit_full         = ifelse(outside == FALSE, TRUE, FALSE)
, fit_red          = TRUE
, lme4_fit_full    = NULL
, spec_name        = spec_name
, outside          = outside
, uncertainty_list = uncertainty_list)
  
}

## If a previously saved response wasn't found, or previous responses are not being used, calculate new responses
if (!exists("saved_spec_res")) {

## recover phylogeny and branch matrices for full and reduced model
z_mats_phylos <- recov_spec_z(
  full_phylo = temp_phy
, species    = spec_name)

## Estimate for the absent species
  if (i == 1) {
    
new_species_estimate_temp <- estiamte_response_unknown_spec(  
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full         
, species           = spec_name
, species_body_size = community_spec_bs
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses)

new_species_estimate <- new_species_estimate_temp

  } else {
    
new_species_estimate_temp <- estiamte_response_unknown_spec(  
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full
, species           = spec_name
, species_body_size = community_spec_bs
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses)
    
new_species_estimate <- rbind(new_species_estimate, new_species_estimate_temp)
    
  }

## If a previously saved response was found, combine species' responses
} else {
  
  if (i == 1) {
    
  new_species_estimate <- saved_spec_res
    
  } else {
  
  new_species_estimate <- rbind(
    new_species_estimate
  , saved_spec_res)
  
  }
  
}

})[3]

## record progress every 20 species
if ((i / 20) %% 1 == 0) {

print(paste(paste("species:", i, ",", spec_name, sep = " ")
  , paste("Round time = ", paste(round(time_check[i], 3), "seconds", sep = " "), sep = " ")
  , paste("Total time = ", paste(round(sum(time_check), 3), "seconds", sep = " "), sep = " ")
      ## mean and sd of the most recently estimated species at day 2
  , paste("Quantiles: day 2 = "
    , paste(as.character(round(quantile(unlist(new_species_estimate[(nrow(new_species_estimate)) - 6, -(1:4)])
    , probs = c(0.10, 0.50, 0.90)), 4))
    , collapse = "_"), sep = " ")
  , sep = " -- "))

}

}

print("out of first loop: unknown respones")

## Temp for debugging 
response_write_logcial(
  which_model          = which_model
, new_species_estimate = new_species_estimate
, prev_responses       = NULL
, uncertainty_list     = uncertainty_list)
  
## also estimate responses for the species with data
for (j in 1:nrow(spec_with_dat)) {
 
  ## once again, check if previous responses have been saved
  if (!is.null(prev_responses)) {
   saved_spec_res <- response_check_logcial(
     prev_responses = prev_responses
  ,  species        = spec_with_dat[j, 1]
   )
   ## If a species response wasn't found, remove this saved variable
   if (nrow(saved_spec_res) == 0) {
     rm(saved_spec_res)
   }
  }
  
  if (!exists("saved_spec_res")) {
   
  if (j == 1) {  

have_species_estimate <- estimate_response_known_spec(  
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full                          ## data.full used for outside == TRUE or FALSE
, species           = spec_with_dat[j, 1]
, species_body_size = community_spec_bs
, which_model       = which_model
, outside           = outside
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses)

  } else {
  
have_species_estimate <- rbind(have_species_estimate
   , estimate_response_known_spec(  
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full
, species           = spec_with_dat[j, 1]
, species_body_size = community_spec_bs
, which_model       = which_model
, outside           = outside
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses))

    }
  } else {
    
    if (j == 1) {
     
    have_species_estimate <- saved_spec_res      
       
    } else {
      
    have_species_estimate <- rbind(
    have_species_estimate
  , saved_spec_res)  
      
    }
  }
  print(spec_with_dat[j, 1])
  }

print("out of second loop: `known' respones")

## Finally, predict for species missing from the phylogeny or missing from the phylogeny and missing body size info
  ## predict for a species missing from the phylogeny, for now just at the mean of the model, no species level random effect
if (nrow(which_species_no_phylo) != 0) {

for (j in 1:length(which_species_no_phylo)) {
  
    if (!is.null(prev_responses)) {
   saved_spec_res <- response_check_logcial(
     prev_responses = prev_responses
  ,  species        = which_species_no_phylo[j]
   )
   ## If a species response wasn't found, remove this saved variable
   if (nrow(saved_spec_res) == 0) {
     rm(saved_spec_res)
   }
  }
  
    if (!exists("saved_spec_res")) {
  
  if (j == 1) {  

phylo_missing_species_estimate <- estiamte_response_missing_spec(
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full
, species           = which_species_no_phylo[j]
, species_body_size = community_spec_bs[community_spec_bs[["species"]] == which_species_no_phylo[j], 3]
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses)

  } else {
  
phylo_missing_species_estimate <- rbind(phylo_missing_species_estimate
   , estiamte_response_missing_spec(
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full
, species           = which_species_no_phylo[j]
, species_body_size = community_spec_bs[community_spec_bs[["species"]] == which_species_no_phylo[j], 3]
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses))
   
  }
      
    } else {
      
  if (j == 1) { 
    
phylo_missing_species_estimate <- saved_spec_res 
    
  } else {
    
phylo_missing_species_estimate <- rbind(
  phylo_missing_species_estimate
, saved_spec_res     
)
    
  }
    }
  
print(paste(paste("species:", j, ",", which_species_no_phylo[j], sep = " ")
  , paste("Est = ", round(phylo_missing_species_estimate[["model_est"]][nrow(phylo_missing_species_estimate)], 5), sep = " ")
  , sep = " -- "))
  
}

  }

print("out of third loop: missing responses")

    if (!is.null(prev_responses)) {
   saved_spec_res <- response_check_logcial(
     prev_responses = prev_responses
  ,  species        = "Missing"
   )
   ## If a species response wasn't found, remove this saved variable
   if (nrow(saved_spec_res) == 0) {
     rm(saved_spec_res)
   }
    }

    if (!exists("saved_spec_res")) {

phylo_and_bs_missing_species_estimate <- estiamte_response_missing_spec(
  phylo_data        = z_mats_phylos
, lme4_fits         = lme4_fits
, data.full         = data.full
, species           = "Missing"
, species_body_size = mean(community_spec_bs[["body_size_s"]])
, which_model       = which_model
, uncertainty_list  = uncertainty_list
, other_responses   = other_responses)

    } else {
  
 phylo_and_bs_missing_species_estimate <- saved_spec_res
      
}

if (exists("phylo_missing_species_estimate")) {
spec_est_data.frame <- rbind(
  new_species_estimate
, have_species_estimate
, phylo_missing_species_estimate
, phylo_and_bs_missing_species_estimate
  )
} else {
spec_est_data.frame <- rbind(
  new_species_estimate
, have_species_estimate
, phylo_and_bs_missing_species_estimate
  ) 
}

## If running with the consensus phylogeny, overwrite saved responses (if desired)
if (write_resp == TRUE & !exists("saved_spec_res")) {
response_write_logcial(
  which_model          = which_model
, new_species_estimate = spec_est_data.frame
, prev_responses       = NULL
, uncertainty_list     = uncertainty_list)
print("Responses exported")
} else if (write_resp == TRUE & exists("saved_spec_res")) {
response_write_logcial(
  which_model          = which_model
, new_species_estimate = spec_est_data.frame
, prev_responses       = prev_responses
, uncertainty_list     = uncertainty_list)
print("Responses exported")  
}

return(spec_est_data.frame)
  
}

##########
### Misc Functions (low level functions)
##########
  ## Recover the desired species vector from the Z matrix
recov_spec_z      <- function (full_phylo, species) {
  #####
  ## (1) Establish reduced phylogeny
  #####
  
try_phy1 <- drop.tip(full_phylo, species, collapse.singles = TRUE)
  
  #####
  ## (2) calculate full and reduced branch lengths and Z matrix 
  #####
  
## Branchlengths and Z matrix for full phylogeny
brl_mat <- compute.brlen(full_phylo)
phyZ    <- phylo.to.Z(brl_mat, stand = FALSE)

## Branchlengths and Z matrix for reduced phylogeny
brl_mat1 <- compute.brlen(try_phy1)
phyZ1    <- phylo.to.Z(brl_mat1, stand = FALSE)
  
  #####
  ## (3) find and store row of interest from full Z matrix by searching for species X
  #####
  
## Row of branch lengths for desired species
spec_est <- phyZ[which(phyZ@Dimnames[[1]] == species),  ]

## Z matrix with removed species
temp_Z   <- phyZ[-which(phyZ@Dimnames[[1]] == species), ]
  
  #####
  ## (4) reduce Z matrix
  #####

## Determine which columns of the matrix have a single value 
## these are the branches that lead to nodes
unique_row_vec <- numeric(ncol(temp_Z))
for (i in 1:ncol(temp_Z)) {
  unique_row_vec[i]   <- ifelse(length(which(temp_Z[, i] > 0)) == 1, 1, 0)
}

## store the sparse matrix as a matrix in order to interact with it more cleanly
temp_Zm               <- as.matrix(temp_Z)
temp_Zm[temp_Zm != 0] <- 1 

## Check if there are duplicated columns in the Z matrix
which.duplicate <- which(duplicated(temp_Zm, MARGIN = 2))

## If there is a duplicated branch, combine the columns of Z (branches in the phylogeny)
  ## Only options for the length are 0 (no duplicates) and 1 (a duplicate)
if (length(which.duplicate) == 1) {
  dup_cols <- which(apply(temp_Zm, 2, all.equal, current = temp_Zm[ , which.duplicate]) == TRUE)
  
## Replace the branch further distant in time with the sum of branch lengths 
last_brl                <- sqrt(rowSums(temp_Z[, dup_cols]^2))
temp_Z[, min(dup_cols)] <- last_brl
last_brl                <- last_brl[which(last_brl != 0)]

## branch lengths leading to the removed species for additional uncertainty
  ## Composed of the first bit of the "duplicated columns" (ancestral piece that is
  ## merged with the most closely related species) + length leading to the removed species
last_brl_uncer <- temp_Z[min(which(temp_Z[, min(dup_cols)] != 0)), min(dup_cols)] + 
  spec_est[max(which(spec_est != 0))]

## Remove the branch most recent in evolutionary time from both Z matrix and species vector
temp_Z    <- temp_Z[, -max(dup_cols)]

## Check for columns that have all 0s
all.zeros <- which(colSums(temp_Z) == 0)

## Remove those columns of all 0s
temp_Z    <- temp_Z[, -all.zeros]

## The absence of a duplicated branch occurs in two scenarios.
 ## 1) The species is an out-group and
 ## 2) The species is part of an unresolved clade (a "star" portion of the phylogeny)
  ## These cases are dealt with in a different way 
} else {

  ## If the species is an out-group it will have only a single branch. In this case the first
   ## two branches must be removed to remove the branch leading to the outgroup and the base
    ## branch to remove the evolutionary change associated with all species from the one outgroup species 
if (length(which(spec_est != 0)) == 1) {
  temp_Z  <- temp_Z[, -c(1, 2)]
} else {

## Check for columns that have all 0s
all.zeros <- which(colSums(temp_Z) == 0)

## Remove those columns of all 0s 
temp_Z    <- temp_Z[, -all.zeros]

  }
  
}

  #####
  ## (5) reduce species vector of Z
  #####

## Need to determine the branches (and thus random effects) that I am ignoring for the species I want to estimate?
  ## and reduce the length of the vector so that it corresponds to the random effect estimates I need

## Currently the random effects for my species include the last branch leading to it *and* the branch that I am
  ## now estimating to its most closely related neighbor (which had two branches collapsed due to my editing
  ## of the Z matrix above (the check for a row that has two columns, where the values in that column are only in that row))

## The question becomes: what is a generic way to remove the value in the column that corresponds to the additional
## column I need to remove that is leading to the most closely related species?
  ## as this will only occur when the group I am estimating for is part of a two species mono-phyletic clade,
    ## when it is an out group this doesn't need to be removed

## This finds the unique branch (most recent in time) to the given species I want to estimate for
  # !! last_b_spec <- as.numeric(names(which(apply(phyZ[, which(spec_est > 0)], 2, function(x) length(which(x > 0))) == 1)))
  
## Remove the appropriate columns (branches) from the species I want to estimate
  ## Should ever only need to remove that last branch leading to the species of interest
 if (length(which.duplicate) == 1) {

   ## Need to remove the column that had a branch lost due to the merging,
    ## which is max(dup_cols)
   ## and the branch that led to the species we just removed, which is the last entry in the 
    ## spec_est vector
   spec_est <- spec_est[-c(max(dup_cols), max(which(spec_est != 0)))]

 } else {

## The absence of a duplicated branch occurs in two scenarios.
 ## 1) The species is an out-group and
 ## 2) The species is part of an unresolved clade (a "star" portion of the phylogeny)
  ## These cases are dealt with in a different way 
   
     ## If the species is an out-group it will have only a single branch
if (length(which(spec_est != 0)) == 1) {
    
    spec_est <- spec_est[-c(1, 2)]
    last_brl_uncer <- brl_mat[["edge.length"]][which(spec_est == 1)]
    
} else {
  
    ## remove the last entry in spec_est and use that as the uncertainty branch
    last_brl_uncer <- spec_est[max(which(spec_est != 0))]
    spec_est       <- spec_est[-max(which(spec_est != 0))]
    
}

 }
  
  #####
  ## (6) return reduced Z and species vector of Z
  #####
    # phyZ = Z matrix for full phylogeny
    # temp_Z = Z matrix for reduced phylogeny
    # spec_est = Z matrix row for desired species to estimate

  return(list(
    full_Z_mat     = phyZ
  , full_phylo     = full_phylo
  , reduced_Z_mat  = temp_Z
  , reduced_phylo  = try_phy1
  , species_vector = spec_est
  , last_brl       = last_brl_uncer))

}
  ## Matches scientific names to tip #s from phylogeny
match_scinames    <- function (needed_species, have_species, phylo) {
  
## tip numbers of the full avian phylogeny that correspond to the species that I need
  needed_tips   <- numeric(length(needed_species))
for (z in seq_along(needed_tips)) {
 needed_tips[z] <- which(phylo[["tip.label"]] == as.character(needed_species[z]))
}

## tip numbers of the full avian phylogeny that correspond to the species that I have data for  
  have_tips   <- numeric(nrow(have_species))
for (z in seq_along(have_tips)) {
 have_tips[z] <- which(phylo[["tip.label"]] == as.character(have_species[, 1])[z])
}
  
## find which species in the data set I do not have data for
  missing_spec <- which(is.na(match(needed_tips, have_tips)))
  
  return(list(
    needed_tips  = needed_tips
  , have_tips    = have_tips
  , missing_spec = missing_spec
  ))
  
}
  ## Species-level competence estimates
comp_summ         <- function (titer_pred, surv_pred, biting_pred, bird_mos_pred
                                 , bird_mos_samps, day, uncertainty_list
                                 , nsamps) {
  
titer_pred  <- titer_pred[order(titer_pred[["species"]]), ]
titer_resp  <- titer_pred[, -c(1:4)]
surv_pred   <- surv_pred[order(surv_pred[["species"]]), ]
surv_resp   <- surv_pred[, -c(1:5)]
biting_pred <- biting_pred[order(biting_pred[["species"]]), ]
biting_resp <- biting_pred[, -c(1:2)]

## If there is no error X_resp will be numeric. Need to be single column matrices
if (class(titer_resp) == "numeric") {
titer_resp  <- matrix(data = titer_resp, ncol = 1, nrow = length(titer_resp))
surv_resp   <- matrix(data = surv_resp, ncol = 1, nrow = length(surv_resp))
biting_resp <- matrix(data = biting_resp, ncol = 1, nrow = length(biting_resp))
}

## dynamic in a really ugly way
host_comp_results <- matrix(nrow = length(unique(titer_pred[["species"]])) *
    length(day) * length(day), ncol = ncol(titer_resp))  

## Corresponds to the number of unique nsamps
  ## When there is no error this is = 1, when there is error this is equal to nsamps
for (i in 1:ncol(titer_resp)) {
  
  host_comp_data.frame <-  titer_pred[, 1:4]
  
  ## Translate titer through the filter of bird to mosquito transition, because of the cap on the importance of titer
    ## And add a survival column
  ## Check if uncertainty in bird to mosquito transmission is to be included
  if (uncertainty_list[["stan_BtoM_model"]] == FALSE) {
    
  ## predict at 26c (note: this is bird to mosquito transmission which only covers a small range of degrees c. 
    ## This is a temp within the range of the data used to fit the model
host_comp_data.frame <- transform(host_comp_data.frame
 , titer        = titer_resp[, i]
 , transmission = plogis(bird_mos_pred[1, 6] * titer_resp[, i] + bird_mos_pred[2, 6] * 26)
 , survival     = surv_resp[, i])

  } else {
    
host_comp_data.frame <- transform(host_comp_data.frame
 , titer        = titer_resp[, i]
 , transmission = plogis(sample(samps_bird_mos[ , 1], 1) * titer_resp[, i] + sample(samps_bird_mos[ , 2], 1) * 26)
 , survival     = surv_resp[, i])

}

## weight transmission by survival
  host_comp_data.frame <- transform(
    host_comp_data.frame
  , realized_transmission   = transmission * survival
  , realized_titer          = titer * survival)
  
## sort species alphabetically in both data frames
host_comp_data.frame <- host_comp_data.frame[order(as.character(host_comp_data.frame[["species"]])), ]
biting_pred          <- biting_pred[order(as.character(biting_pred[["species"]])), ]

check_spec_match <- data.frame(unique(host_comp_data.frame[["species"]]), unique(biting_pred[["species"]]))

if (length(match(check_spec_match[, 1], check_spec_match[, 2])) != length(unique(host_comp_data.frame[["species"]]))) {
  ## gives error about break, but breaks from function at appropriate time with printed statement and error
    ## so basically works as intended
  print("Error: species dont match"); break 
}

## organize the species in both the titer and survial data frame with the biting data frame
matched_rows <- match(host_comp_data.frame[["species"]], biting_pred[["species"]])
host_comp_data.frame <- transform(host_comp_data.frame
  , biting_pref = biting_resp[matched_rows[!is.na(matched_rows)], i] / 1000)
host_comp_data.frame <- transform(host_comp_data.frame
  , weighted_trans = biting_pref * realized_transmission
  , run            = rep(i, nrow(host_comp_data.frame)))

## stitch together results in the most reasonable way I can think of for speed
if (i == 1) {
  host_comp_data.frame      <- melt(host_comp_data.frame, c("day", "log_dose", "species", "model"))
  host_comp_data.frame_info <- host_comp_data.frame[, 1:5]
  host_comp_temp_res        <- as.matrix(host_comp_data.frame[, -c(1:5)])
} else {
  host_comp_data.frame      <- melt(host_comp_data.frame, c("day", "log_dose", "species", "model"))
  host_comp_temp_res        <- as.matrix(host_comp_data.frame[, -c(1:5)])
}

host_comp_results[, i] <- host_comp_temp_res

if ((i / 100) %% 1 == 0) print(i / nsamps)

}

return(cbind(host_comp_data.frame_info, host_comp_results))

}
  ## Community-level competence estimates using community R0: mosquito to mosquito
 ## Updated to incorporate temperature variation at the scale of the bird communities
comp_prop_sum_R0  <- function (host_comp_summary, bird_prop_dat, comm_detect_est, uncertainty_list
                                 , mosquito_survival, m_to_b_trans_samps, county_temp_data, print_prog = TRUE) {
  
  ## R0 calculation instead of community competence calculation
  
  ## Mosquito to host step:
    ## ( relative abundance of species i ) * 
    ## ( average probability of infecting a mosquito )
      ## ( ((biting_pref * prop) / sum(biting_pref * prop)) ) * 
      ## ( total_transmission / 8 )
  
  ## Host to mosquito step: 
    ## ( relative abundance of species i ) * 
    ## ( probability of infection of a host )
      ## ( prop / mean(prop) ) * 
      ## ( 0.5 )
    
  ## Other parameters:
    ## * ( bite rate on birds = bite rate of each mosquito * ratio of mosquitos to birds )
      ## 0.14 * 10 (from Hamer et al. 2009 and Simpson et al. 2011)
  
  ## First determine if there will be one sample or many (with uncertainty present)
nsamps         <- ifelse(no_uncer == TRUE, 1, uncertainty_list[["samps"]])

  ## Remove the few data rows without a county ID
bird_prop_dat <- bird_prop_dat %>% filter(county != "")
   
  ## match detection estimates with each bird species, adding detection est to data frame of ebird counts
match_vec_de   <- match(bird_prop_dat[["phylo_name"]], comm_detect_est[["species"]])

## Also merge the temperature data into the community dataframe
bird_prop_dat  <- transform(bird_prop_dat, month = as.character(month))
bird_prop_dat  <- left_join(bird_prop_dat, county_temp_data, by = c("county", "year", "month"))
  
detect_pred <- comm_detect_est[, -c(1:3)]
  if (class(detect_pred) == "numeric") {
detect_pred <- matrix(data = detect_pred, ncol = 1, nrow = length(detect_pred))  
  }

## Calculate mosquito to bird transmission with uncertainty
temp_day <- seq(1, 50, by = 1)
temp_LD  <- c(5.5)

## No uncertainty in survival currently...
newdat       <- mosquito_survival %>% filter(Longevity_Days <= max(temp_day))
temp_range   <- seq(min(county_temp_data$temp), max(county_temp_data$temp))

## Mosquito to bird transmission for sample i with uncertainty, use median without uncertainty   
for (j in 1:length(temp_range)) {
  
  if (uncertainty_list[["stan_MtoB_model"]] == TRUE) {
    
temp_samps   <- sample(seq(1, nrow(m_to_b_trans_samps), by = 1), nsamps) 

## Sample the posterior of mosquito to bird transmission, then predict mosquito to bird transmission for 
 ## each day, discounted by mosquito survival on that day

m_to_b_trans_temp <- plogis(
matrix(
    rep(m_to_b_trans_samps[temp_samps, 1], each = length(temp_day))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(m_to_b_trans_samps[temp_samps, 2], each = length(temp_day)) * rep(temp_day, length(temp_samps)) 
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(m_to_b_trans_samps[temp_samps, 3], each = length(temp_day)) * rep(temp_LD, length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(m_to_b_trans_samps[temp_samps, 4], each = length(temp_day)) * rep(temp_range[j], length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(m_to_b_trans_samps[temp_samps, 5], each = length(temp_day)) * 
  rep(temp_day, length(temp_samps)) *
  rep(temp_range[j], length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
))

m_to_b_trans_temp <- cbind(
  data.frame(temp = rep(temp_range[j], length(temp_day)))
, sweep(m_to_b_trans_temp, 1, 
  matrix(newdat[newdat$Temperature == temp_range[j], ][["Survival"]])
  , FUN = "*")
)

## Transmission weighted by mosquito survival
if (j == 1) {
m_to_b_trans <- m_to_b_trans_temp  
} else {
m_to_b_trans <- rbind(m_to_b_trans, m_to_b_trans_temp)
}

  } else {

## just a placeholder for length    
temp_samps   <- sample(seq(1, nrow(m_to_b_trans_samps), by = 1), nsamps) 

## Instead of sampling the posterior, just use the median estimate. The rest of the columns (which are now
 ## just duplicates, are thrown away later)
m_to_b_trans_temp <- plogis(
matrix(
    rep(median(m_to_b_trans_samps[, 1]), each = length(temp_day))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 2]), each = length(temp_day)) * rep(temp_day, length(temp_samps)) 
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 3]), each = length(temp_day)) * rep(temp_LD, length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 4]), each = length(temp_day)) * rep(temp_range[j], length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 5]), each = length(temp_day)) * 
  rep(temp_day, length(temp_samps)) *
  rep(temp_range[j], length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
)
     )

m_to_b_trans_temp <- cbind(
  data.frame(temp = rep(temp_range[j], length(temp_day)))
, sweep(m_to_b_trans_temp, 1, 
  matrix(newdat[newdat$Temperature == temp_range[j], ][["Survival"]])
  , FUN = "*")
)

## Transmission weighted by mosquito survival
if (j == 1) {
m_to_b_trans <- m_to_b_trans_temp  
} else {
m_to_b_trans <- rbind(m_to_b_trans, m_to_b_trans_temp)
}
  
  }
}

  ## Filter out the parameters that are needed from host_comp_summary
host_comp_summary_rt <- host_comp_summary %>% filter(outcome == "realized_transmission")
host_comp_summary_rt <- melt(host_comp_summary_rt, c("day", "log_dose", "species", "model", "outcome"))
host_comp_summary_bp <- host_comp_summary %>% filter(outcome == "biting_pref")
host_comp_summary_bp <- melt(host_comp_summary_bp, c("day", "log_dose", "species", "model", "outcome"))
  ## sum transmission over days, already weighted by bite preference
host_comp_summary_rt <- host_comp_summary_rt %>%
  group_by(species, variable) %>%
  summarize(total_transmission = sum(value))
host_comp_summary_bp <- host_comp_summary_bp %>%
  group_by(species, variable) %>%
  summarize(biting_pref = mean(value)) ## all values are the same so this just collapses
## Spread these summarized data back out to match the structure of m_to_b trans, where
 ## samples from the posterior are arranged as separate columns
host_comp_summary_rt <- spread(host_comp_summary_rt, key = "variable", value = "total_transmission")
host_comp_summary_bp <- spread(host_comp_summary_bp, key = "variable", value = "biting_pref")

## set up these data frames to be able to merge with the community composition data
match_vec_rt <- match(bird_prop_dat[["phylo_name"]], host_comp_summary_rt[["species"]])
match_vec_bp <- match(bird_prop_dat[["phylo_name"]], host_comp_summary_bp[["species"]])

######
## The problem here is that combining all of the samps at a single time produces too large
## of an object for my RAM. Instead, will have to loop over samples, but set up to do as little
## as possible in the loop to maximize efficiency
######

for (i in 1:nsamps) {

check_time   <- system.time({
  
 ## Select a new column of detectabilities
bird_prop_dat       <- transform(bird_prop_dat, detectability = detect_pred[match_vec_de, i])
 
 ## Save new bird_prop_dat to edit for this run
bird_prop_dat_e     <- bird_prop_dat
  
 ## scale the ebird counts using the detection distances
bird_prop_dat_e     <- bird_prop_dat_e %>%
   group_by(county, month, year) %>%
   mutate(detectability = 1 / (detectability / max(detectability)))
 
bird_prop_dat_e     <- transform(bird_prop_dat_e, ebird_prior = ebird_prior * detectability)
 
  ## aggregate ebird data
bird_prop_dat_agg   <- bird_prop_dat_e %>% 
   group_by(county, month, year) %>%
   mutate(prop = ebird_prior / sum(ebird_prior))

bird_prop_dat_agg <- as.data.frame(bird_prop_dat_agg)

## Add the bird species competence data
temp_comm <- cbind(bird_prop_dat_agg
  , data.frame(total_transmission = unlist(host_comp_summary_rt[match_vec_rt, 1 + i]))
  , data.frame(biting_pref        = unlist(host_comp_summary_bp[match_vec_bp, 1 + i]))
  )

## Add the mosquito_to_bird transmission probability matched with the temp in the time and place
m_to_b_trans_sing <- m_to_b_trans[, c(1, i + 1)]
names(m_to_b_trans_sing)[2] <- "trans_prob"
m_to_b_trans_sing <- m_to_b_trans_sing %>% 
  group_by(temp) %>% 
  summarize(tot_m_b_trans = sum(trans_prob))

temp_comm <- left_join(temp_comm, m_to_b_trans_sing, by = c("temp"))

## calculate the community competence using relative biting preference
 sub_comm_comp <- temp_comm %>% 
   group_by(county, month, year) %>%
  mutate(
  ## Starting with an infected mosquito, how many of each host type are infected? A total of "1 host" will become infected,
    ## with relative abundance and biting preference scaling where that 1 infection is localized
   mosquito_to_bird = 
    ## relative abundance weighted by biting preference
     ((biting_pref * prop) / sum(biting_pref * prop)) *       ## bird / bird 
    ## Median estimate for NY99 mosquito to bird transmission, incorporating bird survival, infectious period, and 
     ## mosquito to bird infection probability. No error here, because it is not the focus of the paper, just a reasonable
      ## estimate for this arm of the transmission cycle
    tot_m_b_trans *                                            
    ## mosquito biting rate
      0.14                                                   ## bite / infected mosquito * day
    ## alternative extremely crude estimate for average probability of mosquito to bird infection
      # 0.5 * 30   ## (infected host / bite) * day
    ## First part is: (bird / bird) * (infected host / bite) * (bite / infected mosquito * day) * (day) = (infected bird / infected mosquito)
  ) %>%
   mutate(
  bird_to_mosquito = 
    ## vector of infected birds, scaled by relative biting preference
       ## (bird / bird) * (infected bird / infected mosquito)
     ( (biting_pref * mosquito_to_bird  / sum(biting_pref * mosquito_to_bird)) * sum(mosquito_to_bird)) *    
    ## average transmission probability per bite              
     (total_transmission / 8 ) *                              ## infected mosquito / bite
    ## mosquito biting rate
      0.14 *                                                  ## bite / susceptible mosquito * day
    ## number of infective days
      8 *                                                     ## day       
    ## number of mosquitos per infected bird. Ratio of mosquitos to birds is assumed to be 3 for any community
      3                                                      ## susceptible mosquitos / infected bird
    ## Second part produces a vector of infected mosquito / infected mosquito, attributable to each bird (the sum gives R0): 
     ## (bird / bird) * (infected mosquito / bite) * (bite / susceptible mosquito * day) * (day) * (susceptible mosquitos / infected bird) 
      ## = (infected mosquito / infected bird)
  )

  ## R0 is calculated (mosquito to mosquito) as the product of the two vectors above, which weights the new infection of mosquitoes 
   ## (host_to_mosquito) by the expected number of each infected host that would occur from a bite starting with a single infected 
    ## mosquito (mosquito_to_host)
  sub_comm_comp <- sub_comm_comp %>% 
   group_by(county, month, year) %>% 
    summarize(
   R0       = sum(bird_to_mosquito) 
 , num_spec = length(unique(phylo_name))
 , num_ind  = sum(ebird_prior))  
 
sub_comm_comp <- as.data.frame(sub_comm_comp)
 
  if (i == 1) {
## melt results for storage
comm_comp_res_temp <- melt(sub_comm_comp, c("county", "month", "year"))
## store info columns
comm_comp_info     <- comm_comp_res_temp[, 1:4]
## set up storage matrix
comm_comp_res      <- matrix(data = 0, ncol = nsamps, nrow = nrow(comm_comp_res_temp))
  } else {
## melt results for storage
comm_comp_res_temp <- melt(sub_comm_comp, c("county", "month", "year"))    
  }

## store results from run i
comm_comp_res[, i] <- comm_comp_res_temp[, 5]

})
 
if (i == 1) {
check_time_t <- check_time 
} else {
check_time_t <- check_time_t + check_time
}

len_na  <- length(which(is.na(comm_comp_res_temp[, 3])))
len_nan <- length(which(comm_comp_res_temp[, 3] == "NaN"))

if (print_prog == TRUE) {
print(
  paste("Run #", i, "number_na =", len_na, "number_nan =", len_nan, sep = " "
  , "round time =", run_time = round(check_time, 2)[3]
  , "total_time =", total_time = round(check_time_t, 2)[3])
  )
}

}
  
return(cbind(comm_comp_info, comm_comp_res))

}
  ## Calculate species competence as defined by the effect on community R0 when substituted
   ## Calculates median species effect in each community in which it is found
 ## Updated to incorporate temperature variation at the scale of the bird communities
check_species_importance_no_uncer_R0 <- function (host_comp_summary, bird_prop_dat, comm_detect_est
                                                    , nsamps, mosquito_survival, m_to_b_trans_samps, county_temp_data) {
  
  
 ## Remove the few data rows without a county ID
bird_prop_dat <- bird_prop_dat %>% filter(county != "")
   
 ## match detection estimates with each bird species, adding detection est to data frame of ebird counts
match_vec_de   <- match(bird_prop_dat[["phylo_name"]], comm_detect_est[["species"]])

## Median estimates
bird_prop_dat <- transform(bird_prop_dat, detectability = comm_detect_est[match_vec_de, 2])

## scale the ebird counts using the detection distances
 bird_prop_dat <- bird_prop_dat %>%
   group_by(county, month, year) %>%
   mutate(detectability = 1 / (detectability / max(detectability)))
 bird_prop_dat <- transform(bird_prop_dat, ebird_prior = ebird_prior * detectability)
 
## aggregate ebird data
bird_prop_dat_agg <- bird_prop_dat %>% 
   group_by(county, month, year) %>%
   mutate(prop = ebird_prior / sum(ebird_prior))

## match the two data frames
match_vec <- match(bird_prop_dat_agg[["phylo_name"]], host_comp_summary[["species"]])
temp_comm <- data.frame(bird_prop_dat_agg, host_comp_summary[match_vec, c(2, 3)])

## Aggregate to the species level
temp_comm2 <- temp_comm %>% 
   group_by(county, month, year, phylo_name) %>%
  summarize(
    bite_pref          = mean(bite_pref)
  , weighted_sum       = sum(prop)
  , total_transmission = mean(med_comp))

## Just use median mosquito to bird transmission here, uncertainty wont effect answers here
 ## Calculate mosquito to bird transmission with uncertainty
temp_day  <- seq(1, 50, by = 1)
## No uncertainty in survival currently...
newdat       <- mosquito_survival %>% filter(Longevity_Days <= max(temp_day))
temp_range   <- seq(min(county_temp_data$temp), max(county_temp_data$temp))
temp_LD   <- c(5.5)

for (j in 1:length(temp_range)) {

## Instead of sampling the posterior, just use the median estimate. The rest of the columns (which are now
 ## just duplicates, are thrown away later)
  ## not used here, just set up to match when uncertainty is used
temp_samps   <- sample(seq(1, nrow(m_to_b_trans_samps), by = 1), nsamps)  
  
m_to_b_trans_temp <- plogis(
matrix(
    rep(median(m_to_b_trans_samps[, 1]), each = length(temp_day))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 2]), each = length(temp_day)) * rep(temp_day, length(temp_samps)) 
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 3]), each = length(temp_day)) * rep(temp_LD, length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 4]), each = length(temp_day)) * rep(temp_range[j], length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
) +
matrix(
  rep(median(m_to_b_trans_samps[, 5]), each = length(temp_day)) * 
  rep(temp_day, length(temp_samps)) *
  rep(temp_range[j], length(temp_day) * length(temp_samps))
  , ncol = length(temp_samps)
  , nrow = length(temp_day)
)
     )

m_to_b_trans_temp <- cbind(
  data.frame(temp = rep(temp_range[j], length(temp_day)))
, sweep(m_to_b_trans_temp, 1, 
  matrix(newdat[newdat$Temperature == temp_range[j], ][["Survival"]])
  , FUN = "*")
)

## Transmission weighted by mosquito survival
if (j == 1) {
m_to_b_trans <- m_to_b_trans_temp  
} else {
m_to_b_trans <- rbind(m_to_b_trans, m_to_b_trans_temp)
}

}
  
## Add the mosquito_to_bird transmission probability matched with the temp in the time and place
m_to_b_trans_sing <- m_to_b_trans[, c(1, 2)] ## Doesn't matter which column is used because they are all the medians
names(m_to_b_trans_sing)[2] <- "trans_prob"
m_to_b_trans_sing <- m_to_b_trans_sing %>% 
  group_by(temp) %>% 
  summarize(tot_m_b_trans = sum(trans_prob))

## Also merge the temperature data into the community dataframe and then add the mosquito to bird transmission
 ## estimates, matched with temperature
temp_comm2 <- transform(temp_comm2, month = as.character(month))
temp_comm2 <- left_join(temp_comm2, county_temp_data, by = c("county", "year", "month"))
temp_comm2 <- left_join(temp_comm2, m_to_b_trans_sing, by = c("temp"))

## Calculate weighted competence for each species
temp_comm2 <- temp_comm2 %>% 
   group_by(county, month, year) %>%
  mutate(
  ## Starting with an infected mosquito, how many of each host type are infected? A total of "1 host" will become infected,
    ## with relative abundance and biting preference scaling where that 1 infection is localized
   mosquito_to_bird = 
    ## relative abundance weighted by biting preference
     ((bite_pref * weighted_sum) / sum(bite_pref * weighted_sum)) *       ## bird / bird 
    ## Median estimate for NY99 mosquito to bird transmission, incorporating bird survival, infectious period, and 
     ## mosquito to bird infection probability. No error here, because it is not the focus of the paper, just a reasonable
      ## estimate for this arm of the transmission cycle
    tot_m_b_trans *                                            
    ## mosquito biting rate
      0.14                                                   ## bite / infected mosquito * day
    ## alternative extremely crude estimate for average probability of mosquito to bird infection
      # 0.5 * 30   ## (infected host / bite) * day
    ## First part is: (bird / bird) * (infected host / bite) * (bite / infected mosquito * day) * (day) = (infected bird / infected mosquito)
  ) %>%
   mutate(
  bird_to_mosquito = 
    ## vector of infected birds, scaled by relative biting preference
       ## (bird / bird) * (infected bird / infected mosquito)
     ( (bite_pref * mosquito_to_bird  / sum(bite_pref * mosquito_to_bird)) * sum(mosquito_to_bird)) *    
    ## average transmission probability per bite              
     (total_transmission / 8 ) *                              ## infected mosquito / bite
    ## mosquito biting rate
      0.14 *                                                  ## bite / susceptible mosquito * day
    ## number of infective days
      8 *                                                     ## day       
    ## number of mosquitos per infected bird. Ratio of mosquitos to birds is assumed to be 3 for any community
      3                                                    ## susceptible mosquitos / infected bird
    ## Second part produces a vector of infected mosquito / infected mosquito, attributable to each bird (the sum gives R0): 
     ## (bird / bird) * (infected mosquito / bite) * (bite / susceptible mosquito * day) * (day) * (susceptible mosquitos / infected bird) 
      ## = (infected mosquito / infected bird)
  )

## Determine the rank-importance of each species
temp_comm3 <- temp_comm2 %>%
  group_by(county, month, year) %>%
  arrange(desc(bird_to_mosquito))

## Add max number of species
temp_comm4 <- temp_comm3 %>% 
   group_by(county, month, year) %>%
   summarize(rank = n()) %>%
   right_join(., temp_comm3)

## Give each species a rank based on bird-to-mosquito transmission (proportion of R0 attributable to a given bird)
temp_comm5 <- temp_comm4 %>%
   group_by(county, month, year) %>%
   mutate(spec_rank = seq(1, max(rank), by = 1))

## Give each community a unique number to be able to loop over
temp_comm5[["unique_communities"]] <- temp_comm5 %>%
  group_by(county, month, year) %>%
  group_indices

## black numeric to be filled in later
temp_comm5 <- temp_comm5 %>%
   mutate(spec_dilut_amp_substitutive = numeric(n()))

unique_community   <- unique(temp_comm5[["unique_communities"]])
time_check         <- numeric(length(unique_community))
running_spec_count <- numeric(length(unique_community))

## blank columns for summary stats
temp_comm5 <- transform(temp_comm5
  , full_comm_R0   = numeric(nrow(temp_comm5))
  , minus_spec_R0  = numeric(nrow(temp_comm5))
  , over_1_below_1 = numeric(nrow(temp_comm5))
  , below_1_over_1 = numeric(nrow(temp_comm5)))

for (i in 1:length(unique_community)) {

## subset to the community of interest
mini_com <- droplevels(subset(temp_comm5, unique_communities == unique_community[i]))
    
## cumulative sum of species analyzed for time tracker
running_spec_count[i] <- mini_com[["rank"]][1]

## If this combination exists
  if (nrow(mini_com) > 0) {
      
## R0 for this specific community
  mini_com <- transform(mini_com, full_comm_R0 = sum(bird_to_mosquito))  
      
  for (l in 1:nrow(mini_com)) {
        
  ## community with 1 removed species
  temp_sumcum     <- mini_com[-l, ]
  temp_sumcum_out <- mini_com[l, ]
  
  ## competence of new community with a species removed and species adjusted according to their proportions
  if (nrow(temp_sumcum) > 0) {
  temp_sumcum     <- transform(temp_sumcum, weighted_sum = weighted_sum / sum(weighted_sum))
  
  ## recalculate R0 without species l
   temp_sumcum    <- temp_sumcum %>% 
  mutate(
  ## Starting with an infected mosquito, how many of each host type are infected? A total of "1 host" will become infected,
    ## with relative abundance and biting preference scaling where that 1 infection is localized
   mosquito_to_bird = 
    ## relative abundance weighted by biting preference
     ((bite_pref * weighted_sum) / sum(bite_pref * weighted_sum)) *       ## bird / bird 
    ## Median estimate for NY99 mosquito to bird transmission, incorporating bird survival, infectious period, and 
     ## mosquito to bird infection probability. No error here, because it is not the focus of the paper, just a reasonable
      ## estimate for this arm of the transmission cycle
    tot_m_b_trans *                                            
    ## mosquito biting rate
      0.14                                                   ## bite / infected mosquito * day
    ## alternative extremely crude estimate for average probability of mosquito to bird infection
      # 0.5 * 30   ## (infected host / bite) * day
    ## First part is: (bird / bird) * (infected host / bite) * (bite / infected mosquito * day) * (day) = (infected bird / infected mosquito)
  ) %>%
   mutate(
  bird_to_mosquito = 
    ## vector of infected birds, scaled by relative biting preference
       ## (bird / bird) * (infected bird / infected mosquito)
     ( (bite_pref * mosquito_to_bird  / sum(bite_pref * mosquito_to_bird)) * sum(mosquito_to_bird)) *    
    ## average transmission probability per bite              
     (total_transmission / 8 ) *                              ## infected mosquito / bite
    ## mosquito biting rate
      0.14 *                                                  ## bite / susceptible mosquito * day
    ## number of infective days
      8 *                                                     ## day       
    ## number of mosquitos per infected bird. Ratio of mosquitos to birds is assumed to be 10 for any community
      3                                                      ## susceptible mosquitos / infected bird
    ## Second part produces a vector of infected mosquito / infected mosquito, attributable to each bird (the sum gives R0): 
     ## (bird / bird) * (infected mosquito / bite) * (bite / susceptible mosquito * day) * (day) * (susceptible mosquitos / infected bird) 
      ## = (infected mosquito / infected bird)
  )

  ## proportional change in R0 due to removing species X
   mini_com[l, ] <- transform(mini_com[l, ]
     , minus_spec_R0               = sum(temp_sumcum[["bird_to_mosquito"]])
     , spec_dilut_amp_substitutive = sum(temp_sumcum[["bird_to_mosquito"]]) / full_comm_R0)
   
   if (mini_com[l, ][["full_comm_R0"]] > 1) {
   mini_com[l, ] <- transform(mini_com[l, ]
     , over_1_below_1 = ifelse(minus_spec_R0 < 1, 1, 0))
 } else if (mini_com[l, ][["full_comm_R0"]] < 1) {
   mini_com[l, ] <- transform(mini_com[l, ]
     , below_1_over_1 = ifelse(minus_spec_R0 > 1, 1, 0))
 }
    
  } else {
   
   mini_com[l, ] <- transform(mini_com[l, ]
     , minus_spec_R0               = NA
     , spec_dilut_amp_substitutive = NA)
     
  }
      }
 
  ## Fill in the answers back into the main data frame
  temp_comm5[temp_comm5[["unique_communities"]] == unique_community[i], c(13:17)] <- mini_com[, c(13:17)]
      
      }

if ((i / 100) %% 1 == 0) {
  print(
    paste(
    paste((i / length(unique_community))*100, "%", sep = " ")
  , paste("community:", i, sep = " ")
  , paste("Cumulative Species Analyzed = ", round(sum(running_spec_count), 5), sep = " ")
  , sep = " -- "))
}
      
}

return(temp_comm5)
  
}
    ## Calculate accumulation curve of proportion of infected mosuquitos as a function of increasing species richness
## NOT used in the main analysis, just extra summary
calc_species_accum_curve             <- function (species_importance) {

species_importance_sub  <- species_importance %>% filter(Parameter == "prop_mosquitos_infected")
species_importance_sub2 <- species_importance %>% filter(Parameter == "total_transmission")
species_importance_sub3 <- species_importance %>% filter(Parameter == "weighted_sum")
species_importance_sub4 <- species_importance %>% filter(Parameter == "bite_pref")
species_importance_sub  <- melt(species_importance_sub, c("Scientific_Name", "Parameter"))
species_importance_sub2 <- melt(species_importance_sub2, c("Scientific_Name", "Parameter"))
species_importance_sub3 <- melt(species_importance_sub3, c("Scientific_Name", "Parameter"))
species_importance_sub4 <- melt(species_importance_sub4, c("Scientific_Name", "Parameter"))
  
  ## summarize species importance
species_importance_prop_mos <- species_importance_sub %>% 
  group_by(Scientific_Name) %>%
  summarize(med_comp  = quantile(value, probs = 0.5))

species_importance_prop_mos2 <- species_importance_sub2 %>% 
  group_by(Scientific_Name) %>%
  summarize(med_trans  = quantile(value, probs = 0.5))

species_importance_prop_mos3 <- species_importance_sub3 %>% 
  group_by(Scientific_Name) %>%
  summarize(bird_prop = mean(value))

species_importance_prop_mos4 <- species_importance_sub4 %>% 
  group_by(Scientific_Name) %>%
  summarize(bite_pref = mean(value))

species_importance_prop_mos <- cbind(species_importance_prop_mos
  , species_importance_prop_mos2[, -1])

species_importance_prop_mos <- transform(species_importance_prop_mos
  , bird_prop   = species_importance_prop_mos3[["bird_prop"]]
  , bite_pref   = species_importance_prop_mos4[["bite_pref"]])

  ## Sort species by their importance
species_importance_prop_mos <- species_importance_prop_mos[order(species_importance_prop_mos[["med_comp"]]
  , decreasing = TRUE), ]
species_importance_prop_mos[["Scientific_Name"]] <- factor(species_importance_prop_mos[["Scientific_Name"]]
  , levels = unique(species_importance_prop_mos[["Scientific_Name"]]))

## First adjust medians so that they sum to 1
species_importance_prop_mos <- transform(species_importance_prop_mos
  , med_comp = med_comp / sum(med_comp))

## Then calculate cumulative sum
species_importance_prop_mos <- transform(species_importance_prop_mos
  , cum_prop_mosquitos_infected = cumsum(med_comp)
  , rank                        = seq(1, nrow(species_importance_prop_mos), by = 1))

return(species_importance_prop_mos )

}
