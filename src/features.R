# =============================================================================
# src/features.R
# -----------------------------------------------------------------------------
# Assemblage des couches (cible + covariables) sur la grille hexagonale.
#
# C'est la couche "PRÉPARATION DES DONNÉES" partagée par :
#   - l'outil d'évaluation des covariables (src/covariate_evaluation.R),
#   - le pipeline de modélisation (mod_spa_temp.R).
#
# Elle transforme des sources hétérogènes (polygones, rasters, plusieurs années)
# en un TABLEAU RECTANGULAIRE propre : une ligne par hexagone (× année), une
# colonne par variable. Toute l'agrégation spatiale est déléguée à
# aggregate.to.hex() (src/hex_grid.R).
#
# Dépendances : sf, terra, + src/hex_grid.R
# =============================================================================

library(sf)
library(terra)


#' Charger la ZONE D'ÉTUDE depuis un découpage administratif (MRC / municipalité)
#'
#' La zone d'étude sert de gabarit à toute la suite : elle délimite la grille
#' hexagonale ET permet de découper les données lourdes (on ne traite que ce qui
#' tombe dans la zone). Choisir une municipalité plutôt qu'une MRC réduit donc
#' fortement les temps de calcul.
#'
#' Le filtrage se fait DIRECTEMENT à la lecture, via une requête SQL, pour ne PAS
#' charger tout le fichier administratif en mémoire (le SDA du Québec est
#' volumineux). Seules les entités demandées sont lues.
#'
#' @param zone_cfg   Liste issue du YAML (section `study_zone`) :
#'   - source  : chemin du fichier de découpage (ex. "data/SDA.gpkg").
#'   - layer   : couche = niveau administratif. Ex. "mrc_s" (MRC) ou
#'               "munic_s" (municipalités).
#'   - field   : colonne contenant le nom. Ex. "MRS_NM_MRC" (MRC) ou
#'               "MUS_NM_MUN" (municipalités).
#'   - regions : vecteur des noms à retenir (une ou plusieurs entités).
#' @param target_crs CRS projeté cible ; la zone est reprojetée dedans (mètres).
#'
#' @return Un objet sf d'UNE seule géométrie (union des entités retenues), déjà
#'   dans `target_crs`. Renvoie NULL si `zone_cfg` est absent (l'appelant se
#'   rabat alors sur une autre définition de zone, ex. l'emprise de la cible).
load.study.zone <- function(zone_cfg, target_crs = 32198) {
  if (is.null(zone_cfg) || is.null(zone_cfg$source)) return(NULL)

  # Construction de la requête SQL : on double les apostrophes des noms
  # (ex. "L'Islet") pour éviter toute erreur de syntaxe.
  vals <- paste0("'", gsub("'", "''", zone_cfg$regions), "'", collapse = ", ")
  requete <- sprintf("SELECT * FROM %s WHERE %s IN (%s)",
                     zone_cfg$layer, zone_cfg$field, vals)

  z <- sf::st_read(zone_cfg$source, query = requete, quiet = TRUE)
  if (nrow(z) == 0) {
    stop("Zone d'étude vide : vérifier study_zone (layer / field / regions) dans la config.")
  }

  z <- sf::st_zm(z)                      # retire la dimension Z (SDA en 3D)
  z <- sf::st_union(z)                   # fusionne les entités en une seule zone
  z <- sf::st_transform(z, target_crs)   # CRS projeté commun (mètres)
  sf::st_sf(geometry = z)
}


#' Charger une source de covariable décrite dans la config
#'
#' @param cov Liste décrivant UNE covariable, telle qu'écrite dans le YAML :
#'   - type       : "raster" ou "vector"
#'   - path       : chemin du fichier
#'   - field      : (vecteur) nom de la colonne de valeur
#'   - band       : (raster) numéro de bande (défaut 1)
#'   - value_type : "continuous" | "categorical" | "ordinal"
#'   Alternative : `cov$source` peut contenir directement un objet sf/SpatRaster
#'   déjà chargé (cas des données simulées) — dans ce cas `path` est ignoré.
#'
#' @return Un objet sf (vecteur) ou SpatRaster (raster).
load.covariate.source <- function(cov) {
  # Cas données simulées : la source est déjà un objet en mémoire.
  if (!is.null(cov$source)) return(cov$source)

  if (identical(cov$type, "raster")) {
    terra::rast(cov$path)
  } else {
    sf::st_read(cov$path, quiet = TRUE)
  }
}


#' Agréger toutes les covariables sur la grille (variables supposées statiques)
#'
#' Dans ce squelette, les covariables sont traitées comme INVARIANTES dans le
#' temps (une seule valeur par hexagone, valable pour toutes les années). Des
#' covariables variant par année seraient une extension : il suffirait de boucler
#' sur les années comme pour la cible.
#'
#' @param hex          Grille hexagonale (create.hex.grid()).
#' @param covariates   Liste nommée de covariables (format YAML, voir ci-dessus).
#' @param min_coverage Seuil de couverture minimale transmis à aggregate.to.hex().
#'
#' @return Liste avec :
#'   - $hex   : la grille enrichie d'une colonne par covariable.
#'   - $types : vecteur nommé donnant le value_type de chaque covariable
#'              (utile en aval pour choisir le bon test statistique / encodage).
aggregate.covariates <- function(hex, covariates, min_coverage = 0.25) {
  types <- character(0)

  for (nm in names(covariates)) {
    cov <- covariates[[nm]]
    src <- load.covariate.source(cov)
    band <- if (is.null(cov$band)) 1 else cov$band

    hex <- aggregate.to.hex(
      hex, src,
      value_field  = cov$field,
      value_type   = cov$value_type,
      out_name     = nm,
      min_coverage = min_coverage,
      band         = band
    )
    types[nm] <- cov$value_type
    message("  covariable agrégée : ", nm, " (", cov$value_type, ")")
  }

  list(hex = hex, types = types)
}


#' Agréger la cible TBE pour UNE année donnée
#'
#' @param hex          Grille hexagonale.
#' @param target_sf    sf de la cible (tous millésimes confondus).
#' @param target_cfg   Config de la cible : year_field, value_field, value_type,
#'                     ordinal_levels, agg_method (voir YAML).
#' @param year         Année à extraire.
#' @param min_coverage Seuil de couverture.
#'
#' @return Vecteur (une valeur par hexagone) de la cible pour cette année.
#'   Pour une cible ordinale, un facteur ordonné selon `ordinal_levels`.
aggregate.target.year <- function(hex, target_sf, target_cfg, year, min_coverage = 0.25) {
  yf <- target_cfg$year_field
  vf <- target_cfg$value_field

  # Sous-ensemble de l'année demandée.
  sub <- target_sf[st_drop_geometry(target_sf)[[yf]] == year, ]
  if (nrow(sub) == 0) return(rep(NA, nrow(hex)))

  vt <- if (is.null(target_cfg$value_type)) "ordinal" else target_cfg$value_type
  h <- aggregate.to.hex(hex, sub, value_field = vf, value_type = vt,
                        out_name = "__target__", min_coverage = min_coverage)
  vals <- h[["__target__"]]

  # Encodage ordinal si des niveaux sont fournis (garantit l'ordre correct).
  if (!is.null(target_cfg$ordinal_levels)) {
    vals <- factor(vals, levels = target_cfg$ordinal_levels, ordered = TRUE)
  }
  vals
}


#' Construire le tableau hexagone × année (format "panel long")
#'
#' Sortie destinée AU GROUPE 1 (prédiction t+1) : chaque ligne est un hexagone
#' observé une année donnée, avec ses covariables et l'intensité TBE. C'est le
#' format tabulaire attendu par Random Forest, et la base pour construire la
#' cible décalée (t+1) et les variables de voisinage.
#'
#' @param hex        Grille hexagonale.
#' @param target_sf  sf cible.
#' @param target_cfg Config cible.
#' @param covariates Liste de covariables (YAML).
#' @param years      Vecteur des années à inclure.
#' @param min_coverage Seuil de couverture.
#'
#' @return Liste avec :
#'   - $data  : data.frame (hex_id, year, <covariables>, target).
#'   - $hex   : la grille (géométrie) pour les jointures spatiales ultérieures.
#'   - $types : value_type des covariables.
build.hex.panel <- function(hex, target_sf, target_cfg, covariates, years,
                            min_coverage = 0.25) {

  # 1. Covariables (statiques) agrégées une seule fois.
  agg <- aggregate.covariates(hex, covariates, min_coverage)
  hex_cov <- agg$hex
  cov_names <- names(covariates)
  cov_df <- st_drop_geometry(hex_cov)[, c("hex_id", cov_names), drop = FALSE]

  # 2. Cible agrégée pour chaque année, empilée en format long.
  rows <- list()
  for (y in years) {
    tgt <- aggregate.target.year(hex, target_sf, target_cfg, y, min_coverage)
    df <- cov_df
    df$year   <- y
    df$target <- tgt
    rows[[as.character(y)]] <- df
  }
  data <- do.call(rbind, rows)
  rownames(data) <- NULL

  # Les hexagones sans TBE observée cette année-là représentent l'ABSENCE
  # d'infestation : on les code explicitement comme un niveau "Aucun" plutôt
  # que NA, car l'absence est une information (et non une donnée manquante).
  if (is.ordered(data$target) || is.factor(data$target)) {
    lv <- c("Aucun", levels(data$target))
    data$target <- factor(as.character(data$target), levels = lv, ordered = TRUE)
    data$target[is.na(data$target)] <- "Aucun"
  }

  list(data = data, hex = hex_cov, types = agg$types)
}
