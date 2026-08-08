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
# Fleuve St_Laurent
if(!file.exists(cfg$extraction$GRHQ$outpath) | cfg$extraction$GRHQ$force_extraction){
  source("src/data_extraction.R")
  extract.St_Laurent(cfg)
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
remove <- st_read(cfg$study_zone$remove)
study.zone <- st_difference(study.zone, remove)


# Construiction de la grille hexagonale ========================================
source("src/structuration.R")

# Avec les niveau en texte
tbe.grid.chr <- aggregate_tbe_levels_to_hex(study.zone, tbe, cfg)

# Avec les niveaux en entiers
tbe.grid.int <- convert_level(tbe.grid.chr)


# État actuel de l'épidemie ====================================================



# Statistiques historiques =====================================================



# Évolution spatiotemporelle ===================================================
library(geocmeans)
source("src/spattemp_dynamic.R")
source("src/plotting.R")

# LISA
moran.df <- test_W_mat(tbe.grid.int[,"Niveau_2025"])
moran.I.comparison(moran.df)
# Nous en concluons que la matrice Queen est la meilleure pour faire ressortir
# l'autocorrélation spatiale (identique à Rook dans les circonstances)

W.Queen <- poly2nb(tbe.grid.int) |>
  nb2listw(zero.policy = TRUE)

tbe.grid.lisa <- lisa_analysis(tbe.grid.int, W.Queen, cfg)

# Cartographie par année
for (year in cfg$lisa_analysis$start:cfg$lisa_analysis$end) {
  lisa.map.no.legend(tbe.grid.lisa, paste0("typo_", year), year)
}
lisa.map(tbe.grid.lisa, "typo_2025", "2025")


# Geo c-means
tbe.zscore <- st_drop_geometry(tbe.grid.lisa[,names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "z")]])

## Optimisation m
FCM_selection <- select_parameters(
  algo = "FCM",
  data = tbe.zscore,
  k = cfg$geocmeans$k,
  m = seq(1.1, 3, 0.1),
  classidx = TRUE,
  spconsist = FALSE,
  tol = 0.001,
  seed = 456,
  verbose = FALSE
)
m.comparison(FCM_selection)

# Optimisation alpha
SFCM_selection <- select_parameters(
  algo = "SFCM",
  data = tbe.zscore,
  k = cfg$geocmeans$k,
  m = cfg$geocmeans$m,
  nblistw = W.Queen,
  alpha = seq(0, 2, 0.05),
  classidx = TRUE,
  tol = 0.001,
  seed = 456,
  spconsist  = TRUE,
  verbose = FALSE
)
alpha.comparison(SFCM_selection)

# Calcul
SFCM <- SFCMeans(tbe.zscore, W.Queen, 
                 k = cfg$geocmeans$k,
                 m = cfg$geocmeans$m,
                 alpha = cfg$geocmeans$alpha,
                 tol = 0.0001, standardize = FALSE,
                 verbose = FALSE, seed = 456)
calcqualityIndexes(tbe.zscore, SFCM$Belongings, 1.6)

Cartes.SFCM <- mapClusters(tbe.grid.lisa, SFCM$Belongings, undecided = 0.45)
Cartes.SFCM$ClusterPlot
tbe.grid.lisa <- cbind(tbe.grid.lisa, SFCM$Belongings)
names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "X")] <- paste0("groupe_", gsub("\\D", "", names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "X")]))

mean.fields <- names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "Niveau")]
grouping.fields <- names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "groupe")]

groups.df <- compute.group.means(tbe.grid.lisa, grouping.fields, mean.fields, cfg)
group.mean.comparison(groups.df, grouping.fields)


# Modélisation =================================================================



