###########################################################################################
### Calculate how many complete lists are needed to realistically capture the community ###
###########################################################################################

## Run before community_competence.R because some name changing goes on in that script

## Step 1:
 ## For all of the communities for which there are at least 1300 complete lists

## Step 2: 
 ## Subsample over each of these well sampled communities (i) 100 times (j) for X compelte lists (k)

## Step 3: 
 ## Determine what is the error between proportions of each species in the communities
  ## with all lists and with the subset of X lists, as well as the # of missing species

## Step 4: 
 ## Find the point of diminishing returns

## requires the FULL lists with all meta data intact
 ## Can be loaded if previously saved, or used from "data_clean_ebird.R"
# samp_data_com_full  <- read.csv(paste(ebird_file_name, "FULL", sep = "_"), sep = "")

## use communities with 1300 lists to resample from
best_sampled_com    <- effort_metadata_com %>% filter(num_lists > 1300)

## number of communities to sample over (i)
num_communities <- seq(1, nrow(best_sampled_com))

## number of lists (j)
num_lists    <- seq(5, 120, by = 5)

## number of resampling events (k)
num_resamps  <- seq(1, 100)

## storage of results
com_resampling_results <- expand.grid(
  num_communities
, num_lists
, num_resamps)
names(com_resampling_results) <- c("community", "lists", "sampling_event")
com_resampling_results <- transform(com_resampling_results
  , RMSE              = 0
  , prop_missing_spec = 0
  , num_com           = 0)

 ## Subsample over each of these well sampled communities (i) 100 times (j) for X compelte lists (k)
for (i in 1:nrow(best_sampled_com)) {
  
  temp_county <- best_sampled_com[i, ][["county"]]
  temp_month  <- best_sampled_com[i, ][["month"]]
  temp_year   <- best_sampled_com[i, ][["year"]]
  
  ## Extract the appropriate list from samp_data_com
  sing_com <- samp_data_com_full %>% filter(
    COUNTY == temp_county
  , MONTH  == temp_month
  , YEAR   == temp_year)

  ## subsample
   temp_uni <- unique(sing_com[["SAMPLING.EVENT.IDENTIFIER"]])
   
  ## summarize this community to make the comparison to the subsampled community
   sing_com_sum <- sing_com %>% 
     group_by(COUNTY, MONTH, YEAR, SCIENTIFIC.NAME) %>% 
     summarise(total_obs = sum(OBSERVATION.COUNT)) %>% 
     mutate(prop_obs     = total_obs / sum(total_obs))
  
for (j in seq_along(num_lists)) {
  
  ## check how many communities would exist if this were the cutoff
  num_com <- nrow(effort_metadata_com[effort_metadata_com[["num_lists"]] > (num_lists[j] - 1), ])
  
  for (k in seq_along(num_resamps)) {

  ## sample from the unique lists for the current community
  temp_lists <- sample(temp_uni, num_lists[j])
   
  ## subsample the community
  sing_com_temp <- sing_com %>% 
    filter(SAMPLING.EVENT.IDENTIFIER %in% temp_lists)
  
  ## summarize the subset community
  sing_com_temp_sum <- sing_com_temp %>% 
     group_by(COUNTY, MONTH, YEAR, SCIENTIFIC.NAME) %>% 
     summarise(total_obs = sum(OBSERVATION.COUNT)) %>% 
     mutate(prop_obs     = total_obs / sum(total_obs))
  
  ## match species and fill in 0s for missing species
  both_res <- left_join(sing_com_sum[, -c(1, 2, 3)]
    , sing_com_temp_sum[, -c(1, 2, 3)]
    , by = "SCIENTIFIC.NAME")
  both_res[is.na(both_res[["prop_obs.y"]]), ][["prop_obs.y"]] <- 0
  
  ## calculate RMSE and total # of missing species between the two vectors or prop_obs
  resamp_est <- both_res %>% 
    summarise(
      RMSE              = sqrt(sum((prop_obs.x - prop_obs.y)^2))
    , prop_missing_spec = length(which(is.na(total_obs.y))) / length(which(!is.na(total_obs.x))))
  
  ## store results
  com_resampling_results[
    com_resampling_results[["community"]] == i &
    com_resampling_results[["lists"]] == num_lists[j] &
    com_resampling_results[["sampling_event"]] == k
    , c(4, 5, 6)] <- data.frame(resamp_est[1, ], num_com)
  
  }

print(c(i, j))
  
}
  
}

## write RDS for use later
#write_rds(com_resampling_results, "com_resampling_results.rds")

#########
## summarise these results and plot to find the point of diminishing returns
## find the point where there is a good tradeoff between maintaining the largest number of communities
## while still being confident in their ability to represent the bird community
########

com_resampling_results_sum <- com_resampling_results %>% 
  group_by(community, lists) %>%
  summarise(
    num_com     = mean(num_com)
  , RMSE_median = quantile(RMSE, 0.50)
  , RMSE_low    = quantile(RMSE, 0.05)
  , RMSE_high   = quantile(RMSE, 0.95)
  , prop_missing_spec_median = quantile(prop_missing_spec, 0.50)
  , prop_missing_spec_low    = quantile(prop_missing_spec, 0.05)
  , prop_missing_spec_high   = quantile(prop_missing_spec, 0.95))

com_resampling_results_sum_sum <- com_resampling_results %>%
  group_by(lists) %>%
  summarise(
    num_com     = mean(num_com)
  , RMSE_median = quantile(RMSE, 0.50)
  , RMSE_low    = quantile(RMSE, 0.05)
  , RMSE_high   = quantile(RMSE, 0.95)
  , prop_missing_spec_median = quantile(prop_missing_spec, 0.50)
  , prop_missing_spec_low    = quantile(prop_missing_spec, 0.05)
  , prop_missing_spec_high   = quantile(prop_missing_spec, 0.95))

com_resampling_results_sum_sum <- com_resampling_results_sum_sum %>%
  mutate(
    RMSE_diff              = c(0, diff(RMSE_median))
  , prop_missing_spec_diff = c(0, diff(prop_missing_spec_median))
  , prop_diff_com          = c(0, diff(num_com)))

theme_update(axis.text.x = element_text(size = 13))

ggplot(com_resampling_results_sum, aes(as.factor(lists), RMSE_median)) + geom_violin() +
  xlab("Number of Complete Lists") + 
  ylab("Median RMSE for each community")

ggplot(com_resampling_results_sum, aes(lists, RMSE_median)) + geom_point() +
  xlab("Number of Complete Lists") + 
  ylab("Median RMSE for each community")

com_resampling_results_sum <- transform(com_resampling_results_sum,
  num_com = factor(num_com, levels = unique(com_resampling_results_sum$num_com)))

ggplot(com_resampling_results_sum, aes(as.factor(num_com), RMSE_median)) + geom_violin() +
  xlab("Total number of Texas Communities") + 
  ylab("Median RMSE for each community") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(com_resampling_results_sum, aes(log10(num_com), RMSE_median)) + geom_point() +
  xlab("Number of Complete Lists") + 
  ylab("Median RMSE for each community")

ggplot(com_resampling_results_sum, aes(as.factor(lists), prop_missing_spec_median)) + geom_violin() +
  xlab("Number of Complete Lists") + 
  ylab("Proportion of species missing")

ggplot(com_resampling_results_sum, aes(lists, prop_missing_spec_median)) + geom_point() +
  xlab("Number of Complete Lists") + 
  ylab("Proportion of species missing")

ggplot(com_resampling_results_sum, aes(as.factor(num_com), prop_missing_spec_median)) + geom_violin() +
  xlab("Total number of Texas Communities") + 
  ylab("Proportion of species missing") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(com_resampling_results_sum, aes(num_com, prop_missing_spec_median)) + geom_point() +
  xlab("Total number of Texas Communities") + 
  ylab("Proportion of species missing")

## Extra layer of summary
ggplot(com_resampling_results_sum_sum[-1, ], aes(log10(prop_diff_com * - 1), RMSE_diff * - 1)) + 
  geom_point() +
  geom_smooth(se = FALSE) +
  xlab("Gain in number (log10) of Texas Communities") + 
  ylab("Loss in RMSE")

ggplot(com_resampling_results_sum_sum[-1, ], aes(log10(prop_diff_com * - 1), prop_missing_spec_diff * - 1)) + 
  geom_point() +
  geom_smooth(se = FALSE) +
  xlab("Gain in number (log10) of Texas Communities") + 
  ylab("Loss in Species Proportion")
