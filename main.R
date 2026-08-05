library(yaml)
library(sf)

# Ingestion de la config =======================================================
cfg <- read_yaml("config.yaml")


# Extraction des données =======================================================
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


# Lecture des données ==========================================================
tbe <- st_read(cfg$target$source) |>
  st_transform(cfg$project$target_crs)
study.zone <- st_read(cfg$study_zone$source, layer = cfg$study_zone$layer)
study.zone <- subset(study.zone, st_drop_geometry(study.zone)[,cfg$study_zone$field] %in% cfg$study_zone$regions) |>
  st_transform(cfg$project$target_crs)


# Construiction de la grille hexagonale ========================================
source("src/structuration.R")

# Avec les niveau en texte
tbe.grid.chr <- aggregate_tbe_levels_to_hex(study.zone, tbe, cfg)

# Avec les niveaux en entiers
tbe.grid.int <- convert_level(tbe.grid.chr)


# État actuel de l'épidemie ====================================================



# Statistiques historiques =====================================================



# Évolution spatiotemporelle ===================================================
source("src/spattemp_dynamic.R")
source("src/plotting.R")

moran.df <- test_W_mat(tbe.grid.int[,"Niveau_2025"])
moran.I.comparison(moran.df)
# Nous en concluons que la matrice Queen est la meilleure pour faire ressortir
# l'autocorrélation spatiale (identique à Rook dans les circonstances)

W.Queen <- poly2nb(tbe.grid.int) |>
  nb2listw(zero.policy = TRUE)

tbe.grid.lisa <- lisa_analysis(tbe.grid.int, W.Queen, cfg)


# Modélisation =================================================================



