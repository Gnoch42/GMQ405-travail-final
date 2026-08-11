# =============================================================================
# main.R — PROGRAMME PRINCIPAL du projet TBE
# -----------------------------------------------------------------------------
# Déroule TOUTE la méthodologie, dans l'ordre :
#   1. Extraction des données (traitement lourd, fait une seule fois)
#   2. Lecture et structuration (grille hexagonale)
#   3. Ampleur de l'épidémie (cartographie + statistiques annuelles)
#   4. Dynamique spatiale (autocorrélation locale LISA + classification geo-c-means)
#   5. Modélisation (groupe 2 économétrique + groupe 1 prédictif) -> src/modelisation.R
#
# Toutes les sorties (cartes, tableaux, graphiques) sont écrites sous le dossier
# `outputs/`, subdivisé par étape (ampleur/ dynamique/ covariables/ economie/
# prediction/ cartes/). Les paramètres sont centralisés dans config.yaml.
#
# Auteurs : Jean-Yves, Antoine Fortier, Rémy Billette.
# =============================================================================

library(yaml)
library(sf)
library(geocmeans)

# --- Modules du projet ---
source("src/structuration.R")     # grille hexagonale + conversion des niveaux
source("src/spattemp_dynamic.R")  # LISA, matrices de voisinage, moyennes de groupe
source("src/visualisation.R")     # TOUS les graphiques et cartes (style unifié)
source("src/historic_stats.R")    # statistiques d'ampleur annuelles

# Ingestion de la config
cfg <- read_yaml("config.yaml")
out <- cfg$output$root %||% "outputs"   # racine des sorties


# 1. EXTRACTION DES DONNÉES  (traitement lourd, fait une seule fois) --------
# Extraction de la TBE pour les régions d'intérêt.
if (!file.exists(cfg$extraction$TBE$outpath) || cfg$extraction$TBE$force_extraction) {
  source("src/data_extraction.R")
  extract.tbe.data(cfg)
}
# Extraction du fleuve Saint-Laurent (retiré de la zone d'étude).
if (!file.exists(cfg$extraction$GRHQ$outpath) || cfg$extraction$GRHQ$force_extraction) {
  source("src/data_extraction.R")
  extract.St_Laurent(cfg)
}


# 2. LECTURE ET STRUCTURATION DES DONNÉES -----------------------------------
tbe <- st_read(cfg$target$source) |>
  st_transform(cfg$project$target_crs)
# Reconstruction du niveau ordinal à partir de l'indice d'intensité (Ia).
tbe$Niveau <- NA
tbe$Niveau[tbe$Ia == 1] <- "Léger"
tbe$Niveau[tbe$Ia == 2] <- "Modéré"
tbe$Niveau[tbe$Ia == 3] <- "Grave"
tbe$Niveau <- factor(tbe$Niveau, levels = c("Léger", "Modéré", "Grave"), ordered = TRUE)

# Zone d'étude = MRC retenues, moins le fleuve Saint-Laurent.
study.zone <- st_read(cfg$study_zone$source, layer = cfg$study_zone$layer)
study.zone <- subset(study.zone,
                     st_drop_geometry(study.zone)[, cfg$study_zone$field] %in% cfg$study_zone$regions) |>
  st_transform(cfg$project$target_crs)
remove <- st_read(cfg$study_zone$remove)
study.zone <- st_difference(study.zone, remove)

# Grille hexagonale : niveaux en texte puis conversion en entiers (0..3).
tbe.grid.chr <- aggregate_tbe_levels_to_hex(study.zone, tbe, cfg)
tbe.grid.int <- convert_level(tbe.grid.chr)


# 3. AMPLEUR DE L'ÉPIDÉMIE -------------------------------------------------
dir.create(file.path(out, "ampleur"), recursive = TRUE, showWarnings = FALSE)

# Carte de l'intensité observée (année la plus récente) + animation temporelle.
tbe.2014.map <- subset(tbe, ANNEE == 2014) |>
  st_intersection(study.zone) |>
  map.tbe(study.zone)
tmap_save(tbe.2014.map, file.path(out, "ampleur", "tbe_2014.png"))

map.tbe.animation(tbe, study.zone, outdir = file.path(out, "ampleur"))

# Statistiques annuelles (superficie touchée par classe, % du territoire).
historic.stats <- compute.historic.stats(tbe, study.zone)


# 4. DYNAMIQUE SPATIALE ------------------------------------------------------
dir.create(file.path(out, "dynamique"), recursive = TRUE, showWarnings = FALSE)

# --- 4a. Choix de la matrice de pondération spatiale ---
# On compare Rook / Queen / distance inverse via le I de Moran ; la matrice Queen
# maximise l'autocorrélation (identique à Rook sur une grille hexagonale).
moran.df <- test_W_mat(tbe.grid.int[, "Niveau_2025"])
moran.I.comparison(moran.df, export = file.path(out, "dynamique", "moran_matrices.png"))

W.Queen <- poly2nb(tbe.grid.int) |>
  nb2listw(zero.policy = TRUE)

# --- 4b. Typologie LISA (points chauds/froids) par année ---
tbe.grid.lisa <- lisa_analysis(tbe.grid.int, W.Queen, cfg)
for (year in cfg$lisa_analysis$start:cfg$lisa_analysis$end) {
  lisa.map(tbe.grid.lisa, paste0("typo_", year), year, legend = FALSE)
}
lisa.map(tbe.grid.lisa, "typo_2025", "2025", legend = TRUE)

# --- 4c. Classification geo-c-means (typologie temporelle de sévérité) ---
tbe.zscore <- st_drop_geometry(
  tbe.grid.lisa[, names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "z")]])

# Optimisation du paramètre de flou m (inertie + silhouette).
FCM_selection <- select_parameters(
  algo = "FCM", data = tbe.zscore, k = cfg$geocmeans$k,
  m = seq(1.1, 3, 0.1), classidx = TRUE, spconsist = FALSE,
  tol = 0.001, seed = 456, verbose = FALSE)
m.comparison(FCM_selection, export = file.path(out, "dynamique", "geocmeans_m.png"))

# Optimisation du poids spatial alpha (silhouette + incohérence spatiale).
SFCM_selection <- select_parameters(
  algo = "SFCM", data = tbe.zscore, k = cfg$geocmeans$k, m = cfg$geocmeans$m,
  nblistw = W.Queen, alpha = seq(0, 2, 0.05), classidx = TRUE,
  tol = 0.001, seed = 456, spconsist = TRUE, verbose = FALSE)
alpha.comparison(SFCM_selection, export = file.path(out, "dynamique", "geocmeans_alpha.png"))

# Classification finale avec les paramètres retenus.
SFCM <- SFCMeans(tbe.zscore, W.Queen, k = cfg$geocmeans$k, m = cfg$geocmeans$m,
                 alpha = cfg$geocmeans$alpha, tol = 0.0001, standardize = FALSE,
                 verbose = FALSE, seed = 456)
calcqualityIndexes(tbe.zscore, SFCM$Belongings, cfg$geocmeans$m)

Cartes.SFCM <- mapClusters(tbe.grid.lisa, SFCM$Belongings, undecided = 0.45)
Cartes.SFCM$ClusterPlot +
  tm_scalebar(breaks = c(0, 20), text.size = 1, position = c(0.7, 0.25))

# Rattachement des probabilités d'appartenance à la grille + moyennes par groupe.
tbe.grid.lisa <- cbind(tbe.grid.lisa, SFCM$Belongings)
names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "X")] <-
  paste0("groupe_", gsub("\\D", "", names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "X")]))

mean.fields     <- names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "Niveau")]
grouping.fields <- names(tbe.grid.lisa)[startsWith(names(tbe.grid.lisa), "groupe")]
groups.df <- compute.group.means(tbe.grid.lisa, grouping.fields, mean.fields, cfg)
group.mean.comparison(groups.df, grouping.fields,
                      export = file.path(out, "dynamique", "geocmeans_groupes.png"))


# 5. MODÉLISATION  (groupe 2 économétrique + groupe 1 prédictif) ---------
# Nécessite le panel de covariables en cache (src/covariate_download.R puis
# src/covariate_build.R). Écrit sous outputs/economie, outputs/prediction,
# outputs/cartes.
source("src/modelisation.R")
run.modelisation("config.yaml")
