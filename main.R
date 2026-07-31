library(yaml)
library(sf)
library(dplyr)

cfg <- read_yaml("config.yaml")


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

tbe <- st_read(cfg$target$source)
study.zone <- st_read(cfg$study_zone$source, layer = cfg$study_zone$layer)
study.zone <- subset(study.zone, st_drop_geometry(study.zone)[,cfg$study_zone$field] %in% cfg$study_zone$regions)

# Construiction de la grille hexagonale ========================================
source("src/structuration.R")

# Avec les niveau en texte
tbe.grid.chr <- aggregate_tbe_levels_to_hex(study.zone, tbe, cfg)

# Avec les niveaux en entiers
tbe.grid.int <- convert_level(tbe.grid)

