###########################################
### Clean Titer Data from meta-analysis ###
###########################################

titercurves <- read.csv("data/Titer.csv", header = TRUE)

## unique line means infection experiment
titercurves <- transform(titercurves
  , unique_line = factor(rep(seq(1, nrow(titercurves)/8), each = 8))
  , Host_Isolate_Date = Host_Isolate_Date - 1999)

## remove all lines with no titer taken (in spreadsheet each study always extended to day 8)
titercurves <- titercurves %>% filter(!is.na(Titer))

titercurves <- transform(titercurves
  , Titer_SD          = ifelse(Titer_SD == 0, 0.1, Titer_SD)
  , Log_Dose          = log10(titercurves$Virus_Dose)
  , Normalized_Weight = 1 / (Titer_SD / Sample_Size_Death) / max(1 / (Titer_SD / Sample_Size_Death)))

## remove other lineages and rows for which there were no surviving birds
titercurves_reduced <- titercurves %>% 
  filter(Virus_Lineage == "B" | Virus_Lineage == "C")
titercurves_reduced <- titercurves_reduced %>% 
  filter(Sample_Size_Death > 0)

## Censored data. Where graphs reported titer less than 1.7, bump titer to 1.7
titercurves_reduced <- transform(titercurves_reduced
  , censoring = ifelse(Titer <= 1.7, 1, 0))
titercurves_reduced <- titercurves_reduced %>% transform(Titer = ifelse(Titer < 1.7, 1.7, Titer))

titercurves_reduced <- titercurves_reduced %>% 
  transform(
      Species = Scientific_Name
    , exp_day = exp(-Day))

titercurves <- droplevels(titercurves)

## list of the avian species for which I have WNV infection profile data
needed_species <- data.frame(Scientific.Name = unique(titercurves_reduced$Scientific_Name))
