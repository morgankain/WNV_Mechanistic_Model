#########################################################################
## all of the cleaning steps of community_competence.R without running ##
## comp_prop_sum_R0, For use in summary stats cleaned if               ##
## run_spat_temp_aggregated == TRUE                                    ##
#########################################################################

  if (no_uncer == TRUE) {
## Remove excess columns if no uncertainty
comm_comp_summary            <- comm_comp_summary[, 1:5]      
  }
  
## Remove any undefined ("") counties if they stuck around...
comm_comp_summary <- comm_comp_summary %>% filter(county != "")
  
## Reorder factors for plotting
comm_comp_summary[["month"]] <- factor(comm_comp_summary[["month"]], levels = seq(1, 12)) 
comm_comp_summary[["year"]]  <- factor(comm_comp_summary[["year"]], levels = sort(unique(comm_comp_summary[["year"]])))  

## Output of comm_comp_summary is in "wide" format, convert to "long" format for plotting and for stats
names(comm_comp_summary)[4]  <- "outcome"  
comm_comp_summary_m          <- melt(comm_comp_summary, c("county", "month", "year", "outcome"))

## organizing community competence for further analysis and plotting
comm_comp_summary_p <- comm_comp_summary_m %>% 
  filter(outcome == "R0") %>%
  group_by(county, month, year) %>%
  summarize(
            min_comp  = quantile(value, probs = 0.1)
          , med_comp  = quantile(value, probs = 0.5)
          , max_comp  = quantile(value, probs = 0.9)
          , cv_comp   = sd(value) / mean(value)
          , mean_comp = mean(value)
          , var_comp  = var(value))  

comm_comp_summary_p_ns <- comm_comp_summary_m %>% 
  filter(outcome == "num_spec") %>%
  group_by(county, month, year) %>%
  summarize(meannum_spec = mean(value))

## combine county data and ebird metadata
comm_comp_summary_p[["meannum_spec"]] <- comm_comp_summary_p_ns[["meannum_spec"]]

names(effort_metadata_com)[1:3] <- c("county", "month", "year")
effort_metadata_com             <- effort_metadata_com[effort_metadata_com[["county"]] != "", ]

comm_comp_summary_p2 <- merge(
    comm_comp_summary_p
  , effort_metadata_com[, c(1, 2, 3, 6)]
  , by = c("county", "month", "year"), all = TRUE)

county_match         <- match(comm_comp_summary_p2[["county"]], county_data[["County"]])
comm_comp_summary_p2 <- cbind(county_data[county_match, -1], as.data.frame(comm_comp_summary_p2))

comm_comp_summary_p2 <- transform(comm_comp_summary_p2
  , year  = as.numeric(as.character(year)))

## Add temperature to the community competence data frame
comm_comp_summary_p2 <- left_join(comm_comp_summary_p2, county_temp_data, by = c("county", "month", "year"))

comm_comp_summary_p2 <- transform(comm_comp_summary_p2
  , month = as.numeric(as.character(month)))

## Subset data for well observed counties
comm_comp_summary_p2_well_observed <- comm_comp_summary_p2 %>% filter(num_lists >= 80)
