##############################
### Load phylogenetic tree ###
##############################

MyTree <- read.nexus(
  paste(
    "trees/"
  , "consensus"
  , ".nexus"
  , sep = ""))
