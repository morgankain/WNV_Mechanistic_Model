##########################################################################################
### Load and clean results from previously run model for Mosquito to bird transmission ###
##########################################################################################

## mosq_bird_trans = transmission probability from mosquitos to birds at 26 degrees celcius, assuming a dose of 5 log titer units,
  ## over the course of a mosquitos life, incorporating mosquito survival

## If the Stan model has been run previously and the user does not want to rerun it, load data
if (load_stan_mos_bird_res == TRUE & file.exists("saved_output/mosquito_to_bird_transmission.Rds")) {

mos_bird_rds            <- readRDS("saved_output/mosquito_to_bird_transmission.Rds")
mos_bird_model_out_summ <- mos_bird_rds[[1]]
mos_bird_pred           <- mos_bird_rds[[2]]
samps_mos_bird          <- mos_bird_rds[[3]]
rm(mos_bird_rds)

## combine the samps from all chains
samps_mos_bird          <- rbind(samps_mos_bird[,1,], samps_mos_bird[,2,]
  , samps_mos_bird[,3,], samps_mos_bird[,4,])

} else {

########
## Mosquito to bird transmission
########
  
incrate  <- read.csv("data/mos_bird_trans.csv", header = TRUE)

incrate  <- incrate %>% filter(!is.na(Sample_Size))

titer_to_inc <- data.frame(
  transmit = subset(incrate
  , Use_Titer   == "N" & 
    Sample_Size == 25 & 
    Citation    == "Moudy et al 2007")[["Adjusted_Percent_Transmitted_Mean"]]
, titer = subset(incrate
  , Use_Titer == "Y" & 
    Sample_Size == 25 & 
    Citation == "Moudy et al 2007")[["Titer_Level_Mean"]])

## Estimating some data from fits for papers that took titer and not prob
logistic_curve <- function(x, k, b, r){
  k / (1 + b * exp(-r * x))
}

mos_prob <- try(nlxb(transmit ~ k / (1 + b * exp(-r * titer))
  , data = titer_to_inc
  , start = list(k = .75, b = 70, r = .5))
  , silent = TRUE)

mos_plot <- data.frame(titer = seq(0, 10, by = 0.01)
  , transmit = logistic_curve(seq(0, 10, by = 0.01)
    , k = mos_prob[["coefficients"]][1]
    , b = mos_prob[["coefficients"]][2]
    , r = mos_prob[["coefficients"]][3]))

## Bit of an ugly loop to match mean transmission from fitted relationship between transmission and titer
for (i in 1:nrow(incrate)) {
  if (incrate[["Use_Titer"]][i] == "Y") {
    incrate[["Adjusted_Percent_Transmitted_Mean"]][i] <- mos_plot[which(round(mos_plot[["titer"]]
      , digits = 2) == round(incrate[["Titer_Level_Mean"]][i], digits = 2)), 2]
    incrate[["Number_Transmitting_Mean"]][i]          <- round(incrate[["Adjusted_Percent_Transmitted_Mean"]][i] * 
        incrate[["Sample_Size"]][i])
  }
}

## Data associated with mosquito to bird transmission with no NAs and no data from Moudy et al. 2007 that inoculated differently
incrate <- incrate %>% 
  filter(!is.na(Adjusted_Percent_Transmitted_Mean)
    , !is.na(incrate[["Sample_Size"]])
    , Time_Series == "Y"
    , Citation != "Moudy et al 2007")

## Needed to reset factors so that they are a continuous seq of numbers when converted to numeric for Stan model
incrate  <- droplevels(incrate) 

## Separate data into with and without Japanese Encephalitis data
incrateJ <- incrate
incrateY <- incrate %>% filter(Strain != "JEV"); incrateY <- droplevels(incrateY)

########
## Stan model
########

## data
mos_bird.data <- 
  with(incrateJ, 
     list(
      "N"          = nrow(incrateJ)
    , "N_Samp"     = Sample_Size
    , "Samp_Max"   = max(Sample_Size)
    , "N_Trans"    = Number_Transmitting_Mean
    , "Trans_Max"  = max(Number_Transmitting_Mean)
    , "Day"        = Days_Post_Infection
    , "LD"         = Log_Dose
    , "Temp_X_Day" = Temperature_C * Days_Post_Infection
    , "Temp"       = Temperature_C
    , "N_VL"       = length(unique(Virus_Lineage))
    , "IE"         = as.numeric(as.factor(Unique_Line))
    , "N_IE"       = length(unique(Unique_Line))
    , "CIT"        = as.numeric(Citation)
    , "N_CIT"      = length(unique(Citation))
    , "VS"         = as.numeric(Vector_Species)
    , "N_VS"       = length(unique(Vector_Species))))

## Run (very slow) stan model
mos_bird_model_out <- stan(
  file    = "stan/Mosquito_to_Bird_ebird.stan"
, data    = mos_bird.data
, iter    = 8000
, thin    = 2
, warmup  = 2000
, refresh = max(8000/100, 1)
, control = list(max_treedepth = 17, adapt_delta = .94),
  chains  = 4)

## Pleasant way to look at convergence of the model
# launch_shinystan(mos_bird_model_out)

## organize model output
detach("package:tidyr", unload = TRUE)
samps_mos_bird          <- extract(mos_bird_model_out, permuted = FALSE)
library(tidyr)
tidy_mos_bird           <- tidy(mos_bird_model_out)
mos_bird_model_out_summ <- summary(mos_bird_model_out)
mos_bird_trans          <- mos_bird_model_out_summ[["summary"]]

saveRDS(list(mos_bird_model_out_summ, mos_bird_trans, samps_mos_bird)
  , file = "saved_output/mosquito_to_bird_transmission.Rds")

## combine the samps from all chains
samps_mos_bird <- rbind(samps_mos_bird[,1,], samps_mos_bird[,2,], samps_mos_bird[,3,], samps_mos_bird[,4,])

}

########
## Mosquito survival
########
mosqsurv      <- read.csv("data/mosqsurv.csv", header = TRUE)

mosqsurv      <- mosqsurv %>% filter(Citation == "Andreadis et al 2014")

mosq_surv     <- glm(cbind(Surviving_Count
  , Sample_Size - Surviving_Count) ~ Temperature * Longevity_Days
  , family = "binomial"
  , data = mosqsurv)

mosqsurv      <- transform(mosqsurv, surv_fitted = fitted(mosq_surv))

### predict mosquito survival data
new_temp_data <- expand.grid(
    Longevity_Days = seq(1, 120, by = 1)
  , Temperature    = seq(max(county_temp_data$temp), min(county_temp_data$temp)))
new_temp_data <- transform(new_temp_data, Survival = plogis(predict(mosq_surv, newdata = new_temp_data)))
