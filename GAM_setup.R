#####################################################
### Organize data for fitting spatio-temporal GAM ###
#####################################################

## Goal is to fit a spatio-temporal model that has the following spatial predictors:
  # (1) mrf (markov random field) smooth for dominant habitat type to model differences in R0 across habitat 
  # (2) Random effect of county to take care of multiple samples for each county and variation within habitat
  # (3) (X, Y) coordinate pairs to take care of other forms of spatial autocorrelation
  # (4) Human population density

########
## Set up the structure for the mrf smooth, using dominant habitat type
########

if (!exists("hab_layer")) {
  hab_layer <- readShapePoly(county_region_type_path)
}

## match county name to what it is called in the poly object
names(comm_comp_summary_p2)[grep("county", names(comm_comp_summary_p2))] <- "CNTY_NM"

## The next few lines are to determine what is the dominant habitat type in each county
comm_comp_summary_p2_r <- comm_comp_summary_p2 %>%
  dplyr::select(-Population, -min_comp, -max_comp, -cv_comp, -mean_comp, -meannum_spec)

comm_comp_summary_p2_m <- melt(comm_comp_summary_p2_r
  , c("Density", "X", "Y", "CNTY_NM", "month", "year", "med_comp"
    , "var_comp", "num_lists", "temp"))

## retain just the dominant habitat type
comm_comp_summary_p2_m_dh <- 
comm_comp_summary_p2_m %>% 
  group_by(CNTY_NM, month, year) %>%
  mutate(hab_typ = seq(1:n())) %>%
  mutate(dom_hab = max(value)) %>%
    dplyr::filter(dom_hab == value)

## match county name to what it is called in the poly object
names(comm_comp_summary_p2_m_dh)[grep("variable", names(comm_comp_summary_p2_m_dh))] <- "REGIONS"
comm_comp_summary_p2_m_dh[["REGIONS"]] <- as.character(comm_comp_summary_p2_m_dh[["REGIONS"]])

## Annoying renaming needed
for (i in 1:nrow(comm_comp_summary_p2_m_dh)) {
  temp_nam <- strsplit(comm_comp_summary_p2_m_dh[["REGIONS"]][i], "[.]")[[1]]
  ## check where an & was removed
  ifand <- which(temp_nam == "")
  if (length(ifand > 0)) {
  temp_nam[min(ifand)] <- "&"
  temp_nam <- temp_nam[-max(ifand)]
  }

  comm_comp_summary_p2_m_dh[["REGIONS"]][i] <- paste(temp_nam, collapse = " ")
  
}

## Name each of the counties in the data according to their numbers in the poly object
county_match <- match(comm_comp_summary_p2_m_dh[["REGIONS"]], hab_layer@data[["REGIONS"]])
which_count  <- hab_layer@data[["NATRGN"]][unique(county_match)]

## Extract the main coordinate matrices for each hab layer. To counteract the fact that there
 ## may be discontinuous sections of the same habitat, gather all partitions of the habitat together. 
count_dim    <- vector("list", length = length(which_count))
for (i in 1:length(which_count)) {
  rows_needed <- which(hab_layer@data[["NATRGN"]] == i)
    for (j in rows_needed) {
        for (k in 1:length(hab_layer@polygons[[j]])) {
      if (k == 1) {
   temp_dim  <- hab_layer@polygons[[j]]@Polygons[[k]]@coords
   temp_dim  <- rbind(temp_dim, c(NA, NA))
   temp_dimc <- temp_dim
      } else {
   temp_dim1 <- hab_layer@polygons[[j]]@Polygons[[k]]@coords
   temp_dim1 <- rbind(temp_dim1, c(NA, NA))
   temp_dimc <- rbind(temp_dimc, temp_dim1)
      }
        }
      if (j == min(rows_needed)) {
   temp_dimcc <- temp_dimc
      } else {
   temp_dimcc <- rbind(temp_dimcc, temp_dimc)   
      }
    }
  count_dim[[i]] <- temp_dimcc
  names(count_dim)[[i]] <- which_count[i]
}

comm_comp_summary_p2_m_dh <- transform(comm_comp_summary_p2_m_dh
  , NATRGN = as.factor(hab_layer@data[["NATRGN"]][county_match]))

## Name the list appropriately
xt <- list(polys = count_dim)
