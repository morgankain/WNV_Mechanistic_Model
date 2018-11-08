#########################################
### Run the spatio-temporal GAM model ###
#########################################

if (no_uncer == FALSE) {

## Fit the model
comm_comp_spatio_temporal_f <- gamm(
    med_comp ~ 
    ## Spatial fixed
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
  , data    = comm_comp_summary_p2_m_dh)

## Fit the model
comm_comp_spatio_temporal_r <- gamm(
    med_comp ~ 
    ## Spatial fixed
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
  , data    = comm_comp_summary_p2_well_observed_m_dh)

} else {
  
  ## Fit the model
comm_comp_spatio_temporal_f <- gamm(
    med_comp ~ 
    ## Spatial fixed
#   s(X, Y) +
    s(NATRGN, bs = "mrf", xt = xt) + 
    s(log(Density)) +
    ## Temporal fixed
    s(month) + 
    s(year)
    ## Spatial random
  , random = list(CNTY_NM = ~1)
    ## Other
  , method = "REML"
  , data   = comm_comp_summary_p2_m_dh)

## Fit the model
comm_comp_spatio_temporal_r <- gamm(
    med_comp ~ 
    ## Spatial fixed
#   s(X, Y) +
    s(NATRGN, bs = "mrf", xt = xt) + 
    s(log(Density)) +
    ## Temporal fixed
    s(month) + 
    s(year)
    ## Spatial random
  , random = list(CNTY_NM = ~1)
    ## Other
  , method = "REML"
  , data   = comm_comp_summary_p2_well_observed_m_dh)
  
}
