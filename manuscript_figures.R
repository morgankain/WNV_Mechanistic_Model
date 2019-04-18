#######################################################################
### Figures and values used (and maybe not used) for the manuscript ###
#######################################################################

######
### NOTE ###
######
### This script is meant to provide some code for exploring the results, it is not
### perfectly cleaned (more of a dump then other scripts because none of the code
### is *required* for fitting and interpreting results)
### That is, not meant to be run with soruce() *will break* or ctrl+A run.
######

######
### NOTE ###
######
### For example, recreating the Figure 1 boxplot requires run_spat_temp_aggregated to
### be set to TRUE, which increases computation time substantially
######

library(extrafont); library(scales); library(tikzDevice)

########
## Dataset info
########

## check number of communities in each data set
nrow(effort_metadata_com)

nrow(
effort_metadata_com %>% filter(num_lists <= 5)
)

## check species coverage
samp_data_with_meta <- left_join(samp_data_com, effort_metadata_com[, c(1,2,3,6)], by = c("county", "month", "year"))

n_distinct(samp_data_with_meta[["phylo_name"]])

samp_data_with_meta <- samp_data_with_meta %>%
  filter(!is.na(num_lists)) %>%
  filter(num_lists >= 80) 

spec_cover <- samp_data_with_meta %>%
  count(phylo_name) 

spec_cover <- spec_cover %>% 
  mutate(prop_com = n / max(n))

## proportion of communities in which a few specific species are found
spec_cover %>% filter(phylo_name == "Cardinalis_cardinalis")
spec_cover %>% filter(phylo_name == "Aphelocoma_ultramarina")
spec_cover %>% filter(phylo_name == "Cyanocitta_cristata")
spec_cover %>% filter(phylo_name == "Zenaida_macroura")

quantile(spec_cover[["prop_com"]], c(0.05, 0.50, 0.95))

length(which(spec_cover[["n"]] < 26))

########
## Crow bite pref
########

host_comp_summary_c <- host_comp_summary %>% 
  filter(species == "Corvus_brachyrhynchos") %>%
  filter(outcome == "biting_pref")

quantile(host_comp_summary_c[1, -c(1:5)], c(0.05, 0.50, 0.95))

########
## Histogram of R0
########

head(comm_comp_summary_p2)

ggplot(comm_comp_summary_p2_well_observed, aes(med_comp)) + 
  geom_histogram(aes(y=..count..), bins = 100) + 
  geom_vline(xintercept = 1, linetype = "dashed", col = "black", lwd = 1) +
#  stat_density(geom = "line") +
  geom_hline(yintercept = 0, colour = "white", size = 1) +
  xlab("Community R0") +
  ylab("Number of Communities") +
  facet_wrap(~month)

ggplot(comm_comp_summary_p2, aes(med_comp)) + 
  geom_histogram(aes(y=..count..), bins = 100) + 
#  geom_vline(xintercept = median(comm_comp_summary_p2[["med_comp"]])
#    , linetype = "dotted", col = "grey", lwd = 1) +
#    geom_vline(xintercept = median(comm_comp_summary_p2_well_observed[["med_comp"]])
#    , linetype = "dotted", col = "blue", lwd = 1) +
  geom_histogram(data = comm_comp_summary_p2_well_observed, aes(y=..count..)
    , colour = "blue", bins = 100, fill = "blue") +
  geom_vline(xintercept = 1, linetype = "dashed", col = "black", lwd = 1) +
#  stat_density(geom = "line") +
  geom_hline(yintercept = 0, colour = "white", size = 1) +
  xlab("Community R0") +
  ylab("Number of Communities") +
  facet_wrap(~month)

########
## Estimates of R0 by space ***prior to fitting the spatial GAM model***
## Plot both the results from the full model and the model with temperature aggregated over space
## Main text present the model with the reduced dataset and all uncertainty
## Supplement present both full dataset with all uncertainty and reduced dataset with no uncertainty
########
comm_comp_summary_p2_well_observed_m_dh_s <- comm_comp_summary_p2_well_observed_m_dh %>%
  group_by(month, REGIONS) %>%
  summarize(R0 = median(med_comp))

## Lines for each ecoregion that are a bit hard to look at
ggplot(comm_comp_summary_p2_well_observed_m_dh_s[
  comm_comp_summary_p2_well_observed_m_dh_s$REGIONS != "COASTAL SAND PLAIN" & 
  comm_comp_summary_p2_well_observed_m_dh_s$REGIONS != "LLANO UPLIFT", ]
  , aes(as.numeric(month), R0)) + 
  geom_line(aes(colour = REGIONS), lwd = 1) +
geom_point(data = comm_comp_summary_p2_well_observed_m_dh_s[
  comm_comp_summary_p2_well_observed_m_dh_s$REGIONS == "COASTAL SAND PLAIN" | 
  comm_comp_summary_p2_well_observed_m_dh_s$REGIONS == "LLANO UPLIFT", ]
  , aes(as.numeric(month), R0, colour = REGIONS), lwd = 3) + 
scale_x_continuous(breaks = c(2, 5, 8, 10, 12), labels = c("Feb", "May", "Aug", "Oct", "Dec")) +
  geom_vline(xintercept = 10, linetype = "dotted", lwd = 0.5) +
  xlab("Month") + ylab(expression('Median R'[0]))

######
## Make a data frame with full results saved as XXX_s1 and 
## one of the results from the aggregated dataset as XXX_s2
######
comm_comp_summary_p2_well_observed_m_dh_s2 <- comm_comp_summary_p2_well_observed_m_dh
comm_comp_summary_p2_well_observed_m_dh_s2 <- transform(comm_comp_summary_p2_well_observed_m_dh_s2, month = as.factor(month))
comm_comp_summary_p2_well_observed_m_dh_s1 <- transform(comm_comp_summary_p2_well_observed_m_dh_s1, model = "Full")
comm_comp_summary_p2_well_observed_m_dh_s2 <- transform(comm_comp_summary_p2_well_observed_m_dh_s2, model = "Reduced")

comm_comp_summary_p2_well_observed_m_dh_sb <- rbind(
  comm_comp_summary_p2_well_observed_m_dh_s1
, comm_comp_summary_p2_well_observed_m_dh_s2)

## Boxplot that is much easier to look at
ggplot(comm_comp_summary_p2_well_observed_m_dh_sb, aes(month, med_comp, colour = model)) + 
  geom_boxplot(lwd = 0.6) +
 scale_colour_manual(
   values = c("steelblue4", "tomato4")
 , name   = "Model"
 , labels = c("Full", "Aggregated")) +
 scale_x_discrete(breaks = c(2, 5, 8, 10, 12), labels = c("Feb", "May", "Aug", "Oct", "Dec")) +
  xlab("Month") + ylab(expression('Median R'[0]))

## Different boxplot for the supplement that are different colors 
ggplot(comm_comp_summary_p2_well_observed_m_dh_sb, aes(month, med_comp, colour = model)) + 
  geom_boxplot(lwd = 0.6) +
 scale_colour_manual(
   values = c("purple4", "springgreen4")
 , name   = "Model"
 , labels = c("Full eBird 
dataset", "Reduced eBird 
dataset")) +
 scale_x_discrete(breaks = c(2, 5, 8, 10, 12), labels = c("Feb", "May", "Aug", "Oct", "Dec")) +
  xlab("Month") + ylab(expression('Median R'[0]))

comm_comp_summary_p2_well_observed_m_dh_s %>% filter(month == 11) %>%
  group_by(CNTY_NM) %>% summarize(R0 = median(med_comp))
comm_comp_summary_p2_well_observed_m_dh_s2 %>%
  group_by(month) %>% summarize(R0 = median(med_comp))

########
## Range of R0, cv etc. !!! (old stuff from before temperature was added) !!!
########

quantile(comm_comp_summary_p2_well_observed[["med_comp"]], c(0.05, 0.50, 0.95))
quantile(comm_comp_summary_p2_well_observed[["cv_comp"]], c(0.05, 0.50, 0.95))

length(which(comm_comp_summary_p2_well_observed[["med_comp"]] < 1.00)) / nrow(comm_comp_summary_p2_well_observed)

## save one of the data frames for comparison between uncertainty propagated and not
# saveRDS(comm_comp_summary_p2_well_observed, file = "comm_comp_summary_p2_well_observed_no_uncer.rds")
# comm_comp_summary_p2_well_observed_no_uncer <- readRDS("comm_comp_summary_p2_well_observed_no_uncer.rds")
# length(which(comm_comp_summary_p2_well_observed_no_uncer[["med_comp"]] > comm_comp_summary_p2_well_observed[["med_comp"]]))

quantile(comm_comp_summary_p2[["med_comp"]], c(0.05, 0.50, 0.95))
quantile(comm_comp_summary_p2[["cv_comp"]], c(0.05, 0.50, 0.95))

length(which(comm_comp_summary_p2[["med_comp"]] < 1.00)) / nrow(comm_comp_summary_p2)

which_at_mean <- which(comm_comp_summary_p2_well_observed[["cv_comp"]] > 
    (quantile(comm_comp_summary_p2_well_observed[["cv_comp"]], c(0.50)) - 0.0005) & 
comm_comp_summary_p2_well_observed[["cv_comp"]] < 
    (quantile(comm_comp_summary_p2_well_observed[["cv_comp"]], c(0.50)) + 0.0005))

comm_comp_summary_p2_well_observed[which_at_mean[2], ]

check_CI <- comm_comp_summary %>% filter(county == "Collin", month == 1, year == 2017)

quantile(check_CI[1, -c(1, 3)], c(0.05, 0.50, 0.95))

length(which(check_CI[1, -c(1, 3)] < 1.00)) 

which_at_mean <- which(comm_comp_summary_p2[["cv_comp"]] > 
    (quantile(comm_comp_summary_p2[["cv_comp"]], c(0.50)) - 0.00005) & 
comm_comp_summary_p2[["cv_comp"]] < 
    (quantile(comm_comp_summary_p2[["cv_comp"]], c(0.50)) + 0.00005))

comm_comp_summary_p2[which_at_mean[1], ]

check_CI <- comm_comp_summary %>% filter(county == "Brewster", month == 11, year == 2013)

quantile(check_CI[1, -c(1, 3)], c(0.05, 0.50, 0.95))

length(which(check_CI[1, -c(1, 3)] < 1.00)) 

########
## Diluters and Amplifiers
########

species_importance_dat_for_plot <- species_importance_p_r[
    species_importance_p_r[["phylo_name"]] == "Zenaida_macroura" | 
    species_importance_p_r[["phylo_name"]] == "Zenaida_asiatica" |
    species_importance_p_r[["phylo_name"]] == "Cardinalis_cardinalis" |
    species_importance_p_r[["phylo_name"]] == "Aphelocoma_ultramarina" |
    species_importance_p_r[["phylo_name"]] == "Cyanocitta_cristata" |
    species_importance_p_r[["phylo_name"]] == "Carpodacus_mexicanus" |
    species_importance_p_r[["phylo_name"]] == "Cyanocorax_yncas"
  , ]

## with robins (used for a presentation)
species_importance_dat_for_plot <- species_importance_p_r[
    species_importance_p_r[["phylo_name"]] == "Zenaida_macroura" | 
    species_importance_p_r[["phylo_name"]] == "Zenaida_asiatica" |
    species_importance_p_r[["phylo_name"]] == "Cardinalis_cardinalis" |
    species_importance_p_r[["phylo_name"]] == "Aphelocoma_ultramarina" |
    species_importance_p_r[["phylo_name"]] == "Cyanocitta_cristata" |
    species_importance_p_r[["phylo_name"]] == "Cyanocorax_yncas" |
    species_importance_p_r[["phylo_name"]] == "Turdus_migratorius"
  , ]

species_importance_dat_for_plot <- species_importance_dat_for_plot[order(species_importance_dat_for_plot[["med_dilut_effect"]]), ]

species_importance_dat_for_plot <- transform(species_importance_dat_for_plot
  , phylo_name = factor(species_importance_dat_for_plot[["phylo_name"]],
          levels = species_importance_dat_for_plot[["phylo_name"]])
  , dilut_type = c(rep("Diluters", 2), rep("Amplifiers", 5)))

## with robins (used for a presentation)
species_importance_dat_for_plot <- transform(species_importance_dat_for_plot
  , phylo_name = factor(species_importance_dat_for_plot[["phylo_name"]],
          levels = species_importance_dat_for_plot[["phylo_name"]])
  , dilut_type = c(rep("Diluters", 3), rep("Amplifiers", 4)))

ggplot(species_importance_dat_for_plot
  , aes(med_dilut_effect, phylo_name, colour = n_com / max(n_com))) + 
  geom_errorbarh(aes(xmin = min_dilut_effect, xmax = max_dilut_effect), height = 0.2, lwd = 1) +
  geom_point(lwd = 3) + 
  geom_vline(xintercept = 1, linetype = "dashed", col = "black", lwd = 0.5) +
  xlab(expression(paste("Proportional Change in", italic(" R")[0], sep = " "))) +
#  xlab(expression(math{R}))
#  xlab(TeX("$\\ensuremath{\\mathcal R}$")) +
  ylab("") +
  #scale_colour_continuous(trans = "reverse", low = muted("red"), mid = "grey", high = muted("blue"), midpoint = 0.50) + 
  scale_colour_gradient2(low = muted("red"), mid = "grey", high = muted("blue"), midpoint = 0.50) +
  guides(color = guide_legend("Proportion of Communities")) + 
  theme(legend.key.size = unit(.55, "cm")
    , legend.position = c(0.2, 0.8)
    , axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 14)
    , axis.title.x = element_text(size = 12)
    , axis.title.y = element_text(size = 14)) +
  scale_y_discrete(labels = 
      c(
   expression(atop(
          "     Northern Cardinal",
          paste("(", italic("Cardinalis cardinalis"), ")", sep = "")))
  ,
   expression(atop(
          "             Transvolcanic jay",
          paste("(", italic("Aphelocoma ultramarina"), ")", sep = "")))
  , 
   expression(atop(
          "                  Blue jay",
          paste("(", italic("Cyanocitta cristata"), ")", sep = "")))
  ,  
   expression(atop(
          "                   House finch",
          paste("(", italic("Carpodacus mexicanus"), ")", sep = "")))
  ,          
   expression(atop(
          "                Green jay",
          paste("(", italic("Cyanocorax yncas"), ")", sep = "")))
  , 
    expression(atop(
          "White-winged dove",
          paste(  "   (", italic("Zenaida asiatica"), ")", sep = "")))
  ,
    expression(atop(
          "      Mourning dove",
          paste("(", italic("Zenaida macroura"), ")", sep = "")))))

ggsave("spec_imp.pdf", width = 8.5, height = 7)

## with robins (used for a presentation)
ggplot(species_importance_dat_for_plot
  , aes(med_dilut_effect, phylo_name, colour = n_com / max(n_com))) + 
  geom_errorbarh(aes(xmin = min_dilut_effect, xmax = max_dilut_effect), height = 0.2, lwd = 1) +
  geom_point(lwd = 3) + 
  geom_vline(xintercept = 1, linetype = "dashed", col = "black", lwd = 0.5) +
  xlab(expression(paste("Proportional Change in", italic(" R")[0], sep = " "))) +
  ylab("") +
  scale_colour_gradient2(low = muted("red"), mid = "grey", high = muted("blue"), midpoint = 0.50) +
  guides(color = guide_legend("Proportion of Communities")) + 
  theme(legend.key.size = unit(.55, "cm")
    , legend.position = c(0.2, 0.8)
    , axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 14)
    , axis.title.x = element_text(size = 12)
    , axis.title.y = element_text(size = 14)) +
  scale_y_discrete(labels = 
      c(
   expression(atop(
          "     Northern Cardinal",
          paste("(", italic("Cardinalis cardinalis"), ")", sep = "")))
  ,
   expression(atop(
          "                  Blue jay",
          paste("(", italic("Cyanocitta cristata"), ")", sep = "")))
  ,       
        
   expression(atop(
          "             Transvolcanic jay",
          paste("(", italic("Aphelocoma ultramarina"), ")", sep = "")))
  , 
   expression(atop(
          "                Green jay",
          paste("(", italic("Cyanocorax yncas"), ")", sep = "")))
  , 
   expression(atop(
          "              American robin",
          paste(  "   (", italic("Turdus migratorius"), ")", sep = "")))     
  ,
    expression(atop(
          "White-winged dove",
          paste(  "   (", italic("Zenaida asiatica"), ")", sep = "")))
  ,
    expression(atop(
          "      Mourning dove",
          paste("(", italic("Zenaida macroura"), ")", sep = "")))))

########
## Species Richness
########

ggplot(comm_comp_summary_p2_well_observed, aes(meannum_spec, med_comp)) + 
  geom_point(lwd = 1) + 
  geom_smooth(lwd = 1, method = "lm") +
  xlab("Species Richness") + 
  ylab(expression(italic(" R")[0]))

########
## Species Abundance and Physiological Competence
########

ggplot(dilut_assump_summary, aes(abundance, mean_trans)) + 
  geom_point(lwd = 1) + 
  geom_smooth(lwd = 1, method = "lm") +
  scale_x_log10() +
  xlab("log(Species Abundance)") + 
  ylab("Species Physiological Competence")

#######
## Spatio Temporal Model
#######

summary(comm_comp_spatio_temporal_f[[2]])
summary(comm_comp_spatio_temporal_r[[2]])
concurvity(comm_comp_spatio_temporal_f[[2]])
concurvity(comm_comp_spatio_temporal_r[[2]])
layout(1); plot(comm_comp_spatio_temporal_f[[2]], scheme = 2, scale = 0)
gam.check(comm_comp_spatio_temporal_r[[2]])
gamres <- resid(comm_comp_spatio_temporal_r[["lme"]], type = "normalized")
better_check(comm_comp_spatio_temporal_r)

par(mfrow = c(2,2))

### predict R0 from the spatial model
  ## reduced data set
comm_comp_summary_p2_well_observed_m_dh <- transform(comm_comp_summary_p2_well_observed_m_dh
  , pred_vals = predict(comm_comp_spatio_temporal_r[[2]]))

  ## complete data set
comm_comp_summary_p2_m_dh <- transform(comm_comp_summary_p2_m_dh
  , pred_vals = predict(comm_comp_spatio_temporal_f[[2]]))

texas_R0_pred_vals <- comm_comp_summary_p2_well_observed_m_dh[, c(4, 5, 6, 10, 12, 14, 15)]

#####
## Each time month gets changed rerun from here
#####
temp_month <- 10

texas_R0_pred_vals_m_s <- texas_R0_pred_vals %>% 
  group_by(REGIONS, month) %>% 
  summarize(
      'Median R0'                   = median(pred_vals, na.rm = TRUE)
   ,  'CV in R0 among communities'  = sd(pred_vals, na.rm = TRUE) / mean(pred_vals, na.rm = TRUE))

#####
## ggplot of variation in R0 across regions over the seasons
#####
names(texas_R0_pred_vals_m_s)[3] <- "Median_R0"
ggplot(texas_R0_pred_vals_m_s[texas_R0_pred_vals_m_s$REGIONS != "COASTAL SAND PLAIN" & texas_R0_pred_vals_m_s$REGIONS != "LLANO UPLIFT", ]
  , aes(month, Median_R0)) + 
  geom_line(aes(colour = REGIONS), lwd = 1) +
geom_point(data = texas_R0_pred_vals_m_s[texas_R0_pred_vals_m_s$REGIONS == "COASTAL SAND PLAIN" | texas_R0_pred_vals_m_s$REGIONS == "LLANO UPLIFT", ]
  , aes(month, Median_R0, colour = REGIONS), lwd = 3) + 
scale_x_continuous(breaks = c(2, 5, 8, 10, 12), labels = c("Feb", "May", "Aug", "Oct", "Dec")) +
  geom_vline(xintercept = 10, linetype = "dotted", lwd = 0.5) +
  xlab("Month") + ylab("Median R0")

## Also try to aggregate by county
texas_R0_pred_vals_m_s <- texas_R0_pred_vals %>% 
  group_by(CNTY_NM, month) %>% 
  summarize(
      'Median R0'                   = median(pred_vals, na.rm = TRUE)
   ,  'CV in R0 among communities'  = sd(pred_vals, na.rm = TRUE) / mean(pred_vals, na.rm = TRUE))

#####
## ggplot of variation in R0 across regions over the seasons
#####
names(texas_R0_pred_vals_m_s)[3] <- "Median_R0"
ggplot(texas_R0_pred_vals_m_s, aes(month, Median_R0)) + 
  geom_boxplot(aes(group = month), lwd = 1) +
scale_x_continuous(breaks = c(2, 5, 8, 10, 12), labels = c("Feb", "May", "Aug", "Oct", "Dec")) +
  geom_vline(xintercept = 10, linetype = "dotted", lwd = 0.5) +
  xlab("Month") + ylab("Median R0")

## Just pick a single month to plot
texas_R0_pred_vals_m_s <- texas_R0_pred_vals_m_s %>% filter(month == temp_month)

## Store a habitat layer for editing each time month is changed (so not to edit the parent data frame) 
hab_layer1month <- hab_layer

hab_layer1month@data <- left_join(hab_layer1month@data, texas_R0_pred_vals_m_s, by = "REGIONS")

hab_layer1month@polygons[[14]] <- NULL

hab_layer1month@data <- hab_layer1month@data[-14, ]

#qtm(shp = hab_layer1month, fill = c("Median_R0"), fill.palette = "Blues", ncol = 2)

tm1 <- tm_shape(hab_layer1month) + tm_fill(col = "CV in R0 among communities", palette = "Greens") +
  tm_borders(col = "black")
tm2 <- tm_shape(hab_layer1month) + tm_fill(col = "Median R0", palette = "Blues") +
  tm_borders(col = "black")
tmap_arrange(tm2, tm1)

#######
## Spatio Temporal Model predictions
#######

summary(comm_comp_spatio_temporal_r[[2]])
summary(comm_comp_spatio_temporal_f[[2]])

as.data.frame(
comm_comp_summary_p2_m_dh %>%
# filter(meannum_spec == 305) %>%
#  filter(year == 2017) %>%
#  filter(county == "Randall") %>%
#  filter(county == "Cameron") %>%
#  filter(county == "Jefferson") %>%
  group_by(month) %>%
#  group_by(REGIONS) %>%
  summarise(median(pred_vals)))

as.data.frame(
comm_comp_summary_p2_well_observed_m_dh %>%
#  filter(meannum_spec == 305) %>%
#  filter(year == 2000) %>%
#  filter(county == "Randall") %>%
#  filter(county == "Cameron") %>%
#  filter(county == "Jefferson") %>%
#  group_by(REGIONS) %>%
  group_by(month) %>%
#  group_by(year) %>%
  summarise(median(pred_vals)))

#######
## Calculate MSE relative to the estimates from the full model
## Take the full model as truth, and check how close the aggregated data models get to the
## predictions from the full model
#######

## First need to make sure all of the data frames are sorted in the correct way, so that the 
 ## correct rows match
comm_comp_summary_p2_well_observed_full <- comm_comp_summary_p2_well_observed_full[
  with(comm_comp_summary_p2_well_observed_full, order(county, month, year)), ]
comm_comp_summary_p2_well_observed_ts <- comm_comp_summary_p2_well_observed_ts[
  with(comm_comp_summary_p2_well_observed_ts, order(county, month, year)), ]
comm_comp_summary_p2_well_observed_tt <- comm_comp_summary_p2_well_observed_tt[
  with(comm_comp_summary_p2_well_observed_tt, order(county, month, year)), ]

#####
## Try some sort of mean absolute difference over the whole time period, or restricted to the
## months that we care about?
#####
MSE_calc <- data.frame(
   model = c("FULL", "TS", "TT", "BCS", "BCT", "EMPTY") 
 , MSE   = c(
 ## full model to itself
   0
 , mean(abs(comm_comp_summary_p2_well_observed_full[["med_comp"]] - comm_comp_summary_p2_well_observed_ts[["med_comp"]])) 
 , mean(abs(comm_comp_summary_p2_well_observed_full[["med_comp"]] - comm_comp_summary_p2_well_observed_tt[["med_comp"]])) 
 , mean(abs(comm_comp_summary_p2_well_observed_full[["med_comp"]] - comm_comp_summary_p2_well_observed_bcs[["med_comp"]])) 
 , mean(abs(comm_comp_summary_p2_well_observed_full[["med_comp"]] - comm_comp_summary_p2_well_observed_bct[["med_comp"]])) 
 ## full model to a mean-only model 
 , mean(abs(comm_comp_summary_p2_well_observed_full[["med_comp"]] - comm_comp_summary_p2_well_observed_mean_model[["med_comp"]]))
 )
)

MSE_calc <- MSE_calc[rev(order(MSE_calc$MSE)), ]
MSE_calc <- transform(MSE_calc, model = factor(model, levels = rev(MSE_calc$model)))
MSE_calc

#####
## ggplot this result
#####

ggplot(MSE_calc
  , aes(model, MSE)) +
  geom_bar(
# aes(fill = uncer), 
  stat = "identity", position = position_dodge()) +
# scale_fill_manual(values = c("steelblue4", "tomato4")) +
  coord_flip() + 
  ylab("Mean Absolute Error") +
#  guides(fill = guide_legend(title = "Model Uncertainty?")) +
  scale_x_discrete(labels = c(
  "Full Model"
, "Temporally Averaged
Bird Community"
, "Spatially Averaged
Bird Community"
,  "Spatially Averaged
Temperature"
, "Mean Model"
, "Temporally Averaged 
Temperature"
    )) +
  xlab("Model") +
  theme(
    panel.grid.major = element_line(colour = "black", size = 0.1)
  , legend.key.size = unit(.55, "cm")
  )

#####
## restricted to the months that we care about?
#####

submonths <- c(4, 5, 9, 10)

MSE_calc_nu <- data.frame(
   model = c("FULL", "TS", "TT", "BCS", "BCT", "EMPTY") 
 , MSE   = c(
 ## full model to itself
   0
 , mean(abs(
   (comm_comp_summary_p2_well_observed_full %>% filter(month %in% submonths))[["med_comp"]] - 
   (comm_comp_summary_p2_well_observed_ts %>% filter(month %in% submonths))[["med_comp"]]
     )) 
 , mean(abs(
   (comm_comp_summary_p2_well_observed_full %>% filter(month %in% submonths))[["med_comp"]] - 
   (comm_comp_summary_p2_well_observed_tt %>% filter(month %in% submonths))[["med_comp"]]
   )) 
 , mean(abs(
   (comm_comp_summary_p2_well_observed_full %>% filter(month %in% submonths))[["med_comp"]] - 
   (comm_comp_summary_p2_well_observed_bcs %>% filter(month %in% submonths))[["med_comp"]]
   )) 
 , mean(abs(
   (comm_comp_summary_p2_well_observed_full %>% filter(month %in% submonths))[["med_comp"]] - 
   (comm_comp_summary_p2_well_observed_bct %>% filter(month %in% submonths))[["med_comp"]]
   )) 
 ## full model to a mean-only model 
 , mean(abs(
   (comm_comp_summary_p2_well_observed_full %>% filter(month %in% submonths))[["med_comp"]] - 
   (comm_comp_summary_p2_well_observed_mean_model %>% filter(month %in% submonths))[["med_comp"]]
   ))
 )
)

MSE_calc_nu <- MSE_calc_nu[rev(order(MSE_calc_nu$MSE)), ]
MSE_calc_nu <- transform(MSE_calc_nu, model = factor(model, levels = MSE_calc_nu$model))
MSE_calc

#####
## ggplot this result
#####
MSE_calc_nu <- transform(MSE_calc_nu, Uncertainty = "None")
MSE_calc    <- transform(MSE_calc, Uncertainty = "All")
MSE_calc_f <- rbind(MSE_calc, MSE_calc_nu)

ggplot(MSE_calc_f
  , aes(model, MSE)) +
#  geom_bar(stat = "identity", position = position_dodge(), fill = "steelblue4") +
  geom_bar(aes(fill = Uncertainty), stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("steelblue4", "tomato4")) +
  coord_flip() + 
  ylab("Mean Absolute Error") +
#  guides(fill = guide_legend(title = "Model Uncertainty?")) +
  scale_x_discrete(labels = rev(c(
 "Full Model"
, "Temporally Averaged
Bird Community"
, "Spatially Averaged
Bird Community"
,  "Spatially Averaged
Temperature"
, "Temporally Averaged 
Temperature"
, "Mean Model"
    ))) +
  xlab("Model") +
  theme(
    panel.grid.major = element_line(colour = "black", size = 0.1)
  , legend.key.size = unit(.55, "cm")
  )

