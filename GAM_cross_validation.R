###########################################################################
### Run the spatio-temporal GAM model using blocked cross-validation:   ###
### leaving out chunks of spatial data and predicting for a subset of   ###
### the values in the left-out piece.                                   ###
### See Roberts et al. 2017: Ecography 40:923-929 for details           ###
###########################################################################

# debug: RMSE_est <- readRDS("RMSE_est_dat.rds")

uni_count <- as.character(droplevels(unique(comm_comp_summary_p2_well_observed_m_dh[["CNTY_NM"]])))
RMSE_est  <- data.frame(
  county   = numeric(0)
, observed = numeric(0)
, pred_f   = numeric(0)
, pred_r1  = numeric(0)
, pred_r2  = numeric(0)
, pred_r3  = numeric(0))

for (i in seq_along(uni_count)) {
  
  ## Remove a county
  temp_dat <- comm_comp_summary_p2_well_observed_m_dh %>%
    filter(CNTY_NM != uni_count[i])
  
  ## Fit a few spatio temporal models of different complexity
  comm_comp_spatio_temporal_f <- gamm(
    med_comp ~ 
    ## Spatial fixed
    s(X, Y) +
    s(NATRGN, bs = "mrf", xt = xt) + 
    s(log(Density)) +
    ## Temporal fixed
    s(month) + 
    s(year)
    ## Spatial random
  , random = list(CNTY_NM = ~1)
    ## Other
  , weights = 1 / var_comp
  , method  = "REML"
  , data    = temp_dat)
  
  comm_comp_spatio_temporal_r1 <- gamm(
    med_comp ~ 
    ## Spatial fixed
    s(X, Y) +
    s(NATRGN, bs = "mrf", xt = xt) + 
 #   s(log(Density)) +
    ## Temporal fixed
    s(month) + 
    s(year)
    ## Spatial random
  , random = list(CNTY_NM = ~1)
    ## Other
  , weights = 1 / var_comp
  , method  = "REML"
  , data    = temp_dat)
    
  comm_comp_spatio_temporal_r2 <- gamm(
    med_comp ~ 
    ## Spatial fixed
    s(X, Y) +
  #  s(NATRGN, bs = "mrf", xt = xt) + 
    s(log(Density)) +
    ## Temporal fixed
    s(month) + 
    s(year)
    ## Spatial random
  , random = list(CNTY_NM = ~1)
    ## Other
  , weights = 1 / var_comp
  , method  = "REML"
  , data    = temp_dat)
  
  comm_comp_spatio_temporal_r3 <- gamm(
    med_comp ~ 
    ## Spatial fixed
  #  s(X, Y) +
    s(NATRGN, bs = "mrf", xt = xt) + 
    s(log(Density)) +
    ## Temporal fixed
    s(month) + 
    s(year)
    ## Spatial random
  , random = list(CNTY_NM = ~1)
    ## Other
  , weights = 1 / var_comp
  , method  = "REML"
  , data    = temp_dat)
  
  ## Data frame for predictions
  newdat <- comm_comp_summary_p2_well_observed_m_dh %>%
    filter(CNTY_NM == uni_count[i])
  
  ## Predict for the left out county using 
  newdat <- newdat %>%
    mutate(
      predvals_f = 
    predict(
     comm_comp_spatio_temporal_f[[2]]
      , newdata = newdat  
      , re.form = NA)
  ,  predvals_r1 = 
    predict(
     comm_comp_spatio_temporal_r1[[2]]
      , newdata = newdat  
      , re.form = NA)
  ,  predvals_r2 = 
    predict(
     comm_comp_spatio_temporal_r2[[2]]
      , newdata = newdat  
      , re.form = NA)
  ,  predvals_r3 = 
    predict(
     comm_comp_spatio_temporal_r3[[2]]
      , newdata = newdat  
      , re.form = NA)
      )
  
  RMSE_est <- rbind(RMSE_est
      , data.frame(
        county   = newdat[["CNTY_NM"]]
      , observed = newdat[["med_comp"]]
      , pred_f   = newdat[["predvals_f"]]
      , pred_r1  = newdat[["predvals_r1"]]
      , pred_r2  = newdat[["predvals_r2"]]
      , pred_r3  = newdat[["predvals_r3"]]))
  
}

RMSE_est   <- droplevels(RMSE_est)

# saveRDS(RMSE_est, file = "RMSE_est_dat.rds")

RMSE_est_c <- RMSE_est %>%
  summarize(
    rmse_f  = sqrt(mean((pred_f - observed)^2))
  , rmse_r1 = sqrt(mean((pred_r1 - observed)^2))
  , rmse_r2 = sqrt(mean((pred_r2 - observed)^2))
  , rmse_r3 = sqrt(mean((pred_r3 - observed)^2)))
  