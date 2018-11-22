##############################################
### Clean Survival Data from meta-analysis ###
##############################################

survival <- read.csv("data/Survival.csv", header = TRUE)

survival <- transform(survival, unique_line = factor(rep(seq(1, nrow(survival)/8), each = 8))
  , Host_Isolate_Date = Host_Isolate_Date - 1999)

## remove the Reisen et al. 2005 that didn't have survival info
survival <- survival %>% filter(!is.na(Titer), !is.na(Survival))

survival <- transform(survival
  , Titer_SD = ifelse(Titer_SD == 0, 0.1, Titer_SD)
  , Log_Dose = log10(survival$Virus_Dose)
  , Normalized_Weight = 1 / (Titer_SD / Sample_Size_Death) / max(1 / (Titer_SD / Sample_Size_Death)))

survival <- droplevels(survival)

## remove other lineages and rows for which there were no surviving birds
survival_reduced <- survival %>% 
  filter(Virus_Lineage == "B" | Virus_Lineage == "C")

survival_reduced <- survival_reduced %>% 
  transform(
      Species = Scientific_Name)

## Clean up survival reduced to be in the right setup for a survival analysis, e.g. daily hazard
survival_reduced <- survival_reduced %>% 
  group_by(unique_line) %>%
  mutate(
    exp_day = exp(-Day)
  , Died = c(0, diff(Alive)) * - 1)
