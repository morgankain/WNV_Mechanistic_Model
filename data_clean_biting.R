################################################
### Clean Biting Data from Hamer et al. 2009 ###
################################################

## Data from Hamer et al. 2009 + eBird prior
mosq_bite_with_prior <- read.delim("data/mosq_bite_with_prior.csv", header = TRUE) 

## sort by proportion
mosq_bite_prop_trial <- mosq_bite_with_prior[order(mosq_bite_with_prior$Proportion), ]

## all non NA values for proportions
spec_prop_comp <- mosq_bite_prop_trial %>% filter(!is.na(Proportion)) %>% dplyr::select(Proportion)

## host species
spec_bites <- mosq_bite_prop_trial %>% filter(!is.na(Bites)) %>% dplyr::select(Bites)

## all non NA values for bites
host_species <- as.numeric(mosq_bite_prop_trial$Scientific_Name)

## proportions must be given as counts for the multinomial, but the total # of observations is not known.
  ## in lieu of this data point, multiply by 1/rarest bird for total counts of each bird
spec_prop_comp <- spec_prop_comp %>% transform(Proportion = round(Proportion / min(filter(spec_prop_comp, Proportion != 0))))

## Prior for bird observations
prev_obs_p <-  mosq_bite_prop_trial %>% dplyr::select(Ebird_Prior)

## Shape and Scale parameter for the prior for biting 
prev_obs_b <- 4

