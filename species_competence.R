####################################
### Calculate Species Competence ###
####################################

## no inclusion of species proportions (species level competence estimates)
host_comp_summary <- comp_summ(
  titer_pred       = comm_titer_est
, surv_pred        = comm_survival_est
, biting_pred      = comm_bite_est
, bird_mos_pred    = bird_mos_pred
, bird_mos_samps   = samps_bird_mos
, day              = seq(1, 8, by = 1)
, uncertainty_list = uncertainty_list
, nsamps           = ifelse(no_uncer == FALSE, uncertainty_list[["samps"]], 1))

  if (no_uncer == TRUE) {
## Remove excess columns if no uncertainty
host_comp_summary        <- host_comp_summary[, 1:6]
  }

## melt in order to summarise in various ways later
host_comp_summary           <- host_comp_summary[host_comp_summary[["variable"]] != "run", ]
names(host_comp_summary)[5] <- "outcome"
host_comp_summary_m         <- melt(host_comp_summary, c("day", "log_dose", "species", "model", "outcome"))

## organizing species competence for plots
host_comp_summary_p <- host_comp_summary_m %>% 
  filter(outcome == "realized_transmission") %>%
  group_by(species, variable) %>%
  summarize(tot_trans = sum(value)); rm(host_comp_summary_m)

## Quantiles
host_comp_summary_p2 <- host_comp_summary_p %>% 
  group_by(species) %>%
  summarize(var_comp = var(tot_trans)
          , min_comp = quantile(tot_trans, probs = 0.1)
          , med_comp = quantile(tot_trans, probs = 0.5)
          , max_comp = quantile(tot_trans, probs = 0.9))

## sort
host_comp_summary_p2 <- host_comp_summary_p2[order(host_comp_summary_p2[["med_comp"]]), ]
host_comp_summary_p2 <- transform( host_comp_summary_p2, species = factor(species, levels = species))

## plot top 30 species
gg_spec_comp <- ggplot(tail(host_comp_summary_p2, 30), aes(species, med_comp)) + 
  geom_point(lwd = 3) +
  geom_errorbar(aes(ymin = min_comp, ymax = max_comp)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
