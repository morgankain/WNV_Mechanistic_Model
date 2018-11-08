#############################################
### Figures and values for the manuscript ###
#############################################

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

ggplot(comm_comp_summary_p2, aes(med_comp)) + 
  geom_histogram(aes(y=..count..), bins = 100) + 
  geom_vline(xintercept = median(comm_comp_summary_p2[["med_comp"]])
    , linetype = "dotted", col = "grey", lwd = 1) +
    geom_vline(xintercept = median(comm_comp_summary_p2_well_observed[["med_comp"]])
    , linetype = "dotted", col = "blue", lwd = 1) +
  geom_histogram(data = comm_comp_summary_p2_well_observed, aes(y=..count..)
    , colour = "blue", bins = 100, fill = "blue") +
  geom_vline(xintercept = 1, linetype = "dashed", col = "black", lwd = 1) +
#  stat_density(geom = "line") +
  geom_hline(yintercept = 0, colour = "white", size = 1) +
  xlab("Community R0") +
  ylab("Number of Communities")

########
## Range of R0, cv etc. (Community R0 paragraph 1)
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
    (quantile(comm_comp_summary_p2_well_observed[["cv_comp"]], c(0.50)) - 0.00005) & 
comm_comp_summary_p2_well_observed[["cv_comp"]] < 
    (quantile(comm_comp_summary_p2_well_observed[["cv_comp"]], c(0.50)) + 0.00005))

comm_comp_summary_p2_well_observed[which_at_mean[1], ]

check_CI <- comm_comp_summary %>% filter(county == "Fort Bend", month == 1, year == 2013)

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
  , dilut_type = c(rep("Diluters", 2), rep("Amplifiers", 4)))

## with robins (used for a presentation)
species_importance_dat_for_plot <- transform(species_importance_dat_for_plot
  , phylo_name = factor(species_importance_dat_for_plot[["phylo_name"]],
          levels = species_importance_dat_for_plot[["phylo_name"]])
  , dilut_type = c(rep("Diluters", 3), rep("Amplifiers", 4)))

library(extrafont); library(scales); library(tikzDevice)

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
layout(1); plot(comm_comp_spatio_temporal_r[[2]], scheme = 2, scale = 0)
gam.check(comm_comp_spatio_temporal_r[[2]])
gamres <- resid(comm_comp_spatio_temporal_r[["lme"]], type = "normalized")
better_check(comm_comp_spatio_temporal_r)

par(mfrow = c(2,2))

texas_R0_pred_vals <- comm_comp_summary_p2_well_observed_m_dh[, c(5, 6, 10, 12, 14, 15)]

texas_R0_pred_vals_m_s <- texas_R0_pred_vals %>% 
  group_by(REGIONS) %>% 
  summarize(
      'Median R0'                   = median(pred_vals)
   ,  'SD in R0 among communities'  = sd(pred_vals))

hab_layer@data <- left_join(hab_layer@data, texas_R0_pred_vals_m_s, by = "REGIONS")

hab_layer1 <- hab_layer
hab_layer1@polygons[[14]] <- NULL

hab_layer1@data <- hab_layer1@data[-14, ]

qtm(shp = hab_layer1, fill = c("Median R0"), fill.palette = "Blues", ncol = 2)

tm1 <- tm_shape(hab_layer1) + tm_fill(col = "SD in R0 among communities", palette = "Greys") +
  tm_borders(col = "black")
tm2 <- tm_shape(hab_layer1) + tm_fill(col = "Median R0", palette = "Blues") +
  tm_borders(col = "black")
tmap_arrange(tm2, tm1)

#######
## Spatio Temporal Model predictions
#######

### predict R0 from the spatial model
  ## reduced data set
comm_comp_summary_p2_well_observed_m_dh <- transform(comm_comp_summary_p2_well_observed_m_dh
  , pred_vals = predict(comm_comp_spatio_temporal_r[[2]]))
  ## complete data set
comm_comp_summary_p2_m_dh <- transform(comm_comp_summary_p2_m_dh
  , pred_vals = predict(comm_comp_spatio_temporal_f[[2]]))

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
  summarise(median(pred_vals))
  )

as.data.frame(
comm_comp_summary_p2_well_observed_m_dh %>%
#  filter(meannum_spec == 305) %>%
#  filter(year == 2000) %>%
#  filter(county == "Randall") %>%
#  filter(county == "Cameron") %>%
#  filter(county == "Jefferson") %>%
  group_by(REGIONS) %>%
#  group_by(month) %>%
#  group_by(year) %>%
  summarise(median(pred_vals))
  )
