library(yaml)

cfg <- read_yaml("ignore/config.yaml")


# Extraction de la donnée pour les régions d'intérêt. Traitement lourd fait une seule fois
if(!file.exists(cfg$extraction$TBE$outpath) | cfg$extraction$TBE$force_extraction){
  source("src/data_extraction.R")
  extract.tbe.data(cfg)
}
# Ambition abandonnée :( Trop lourd
# if(!file.exists(cfg$extraction$IEQM$outpath) | cfg$extraction$IEQM$force_extraction){
#   source("src/data_extraction.R")
#   extract.ieqm.data(cfg)
# }
