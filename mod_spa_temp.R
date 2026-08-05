#########################################################
# mod_spa_temp.R — Pipeline de modélisation spatio-temporelle TBE ----
#########################################################
#
# Programme PRINCIPAL. Il orchestre deux SOUS-PIPELINES distincts (comparés
# séparément, jamais entre eux) :
#
#   GROUPE 1 — PRÉDICTION t+1 (où sera la TBE l'an prochain)
#       Markov spatial · Random Forest/XGBoost · ConvLSTM
#       -> métriques : kappa pondéré, accuracy, RMSE sur rangs
#
#   GROUPE 2 — ANALYSE à un instant t (comprendre les relations)
#       MCO · SLX · SAR · SEM · Durbin spatial · Durbin erreur · GAM
#       -> métriques : AIC/BIC, pseudo-R², test de Moran (diagnostic spatial)
#
# PRÊT À L'EMPLOI : sans configuration, le pipeline tourne sur des DONNÉES
# SIMULÉES (voir config `simulate$use_simulated: TRUE`). Pour passer aux vraies
# données, il suffit de renseigner les chemins dans config.yaml — AUCUNE
# modification de code n'est nécessaire. Les covariables réelles s'ajoutent
# uniquement via la section `covariates:` du YAML.
#
# Lancement :  Rscript mod_spa_temp.R   (ou source() dans RStudio)
#########################################################

rm(list = ls())
graphics.off()

library(sf)
library(yaml)

# --- Modules du projet ---
source("src/hex_grid.R")            # grille + agrégation spatiale générique
source("src/features.R")            # assemblage cible + covariables sur la grille
source("src/models_group1.R")       # prédiction t+1
source("src/models_group2.R")       # analyse à un instant t
# ConvLSTM (torch) : chargé si disponible ; sinon le modèle dégrade proprement.
try(source("src/model_convlstm.R"), silent = TRUE)


# =============================================================================
# CHARGEMENT DES DONNÉES (réelles OU simulées)
# =============================================================================

#' Charger les données du projet selon la config
#'
#' @return Liste : $zone (emprise), $target (sf cible), $covariates (liste),
#'   $target_cfg (config de la cible).
load.project.data <- function(cfg) {

  if (isTRUE(cfg$simulate$use_simulated)) {
    # --- Mode DÉMO : données factices ---
    message(">> Mode DONNÉES SIMULÉES (config simulate$use_simulated = TRUE)")
    source("src/simulate_data.R")
    sim <- simulate.tbe.dataset(
      seed      = cfg$simulate$seed      %||% 42,
      n_years   = cfg$simulate$n_years   %||% 12,
      extent_km = cfg$simulate$extent_km %||% 200
    )
    target_cfg <- list(year_field = "year", value_field = "intensity",
                       value_type = "ordinal",
                       ordinal_levels = c("Léger", "Modéré", "Grave"))
    return(list(zone = sim$zone, target = sim$tbe, covariates = sim$covariates,
                target_cfg = target_cfg))
  }

  # --- Mode RÉEL : lecture depuis les chemins du YAML ---
  message(">> Mode DONNÉES RÉELLES (chemins lus depuis config.yaml)")
  target_crs <- cfg$project$target_crs %||% 32198

  target <- st_read(cfg$target$source, quiet = TRUE)
  # Normalisation des libellés d'intensité (accents/casse incohérents dans la
  # source : "Leger"/"Léger", "Modere"/"Modéré"). On harmonise pour l'ordinal.
  vf <- cfg$target$value_field
  target[[vf]] <- normalize.intensity(target[[vf]], cfg$target$ordinal_levels)
  target <- st_transform(target, target_crs)   # CRS projeté commun (mètres)

  # ZONE D'ÉTUDE : de préférence un découpage administratif (MRC/municipalité)
  # défini dans config `study_zone`. Cela délimite la grille ET allège les
  # traitements (on ne garde que les données de la zone). À défaut, on se rabat
  # sur l'emprise (bounding box) de la cible.
  zone <- load.study.zone(cfg$study_zone, target_crs)
  if (is.null(zone)) {
    message("  (study_zone non défini -> zone = emprise de la cible)")
    zone <- st_as_sfc(st_bbox(target)) |> st_sf(geometry = _)
  } else {
    # Découpage de la cible à la zone : on ne conserve que les entités TBE qui
    # intersectent la zone d'étude, ce qui réduit fortement le volume à traiter.
    target <- target[zone, ]
    message("  Cible découpée à la zone d'étude : ", nrow(target), " entités TBE conservées")
  }

  list(zone = zone, target = target, covariates = cfg$covariates,
       target_cfg = cfg$target)
}


# normalize.intensity() est défini dans src/features.R (partagé avec le build).

# Opérateur "valeur par défaut si NULL" (pratique pour lire le YAML).
`%||%` <- function(a, b) if (is.null(a)) b else a


# =============================================================================
# PRÉPARATION DES DONNÉES DE MODÉLISATION
# =============================================================================

#' Obtenir le panel hex×année (cible + covariables) pour la modélisation
#'
#' Deux sources possibles :
#'   - CACHE de covariables réelles (data/covariates/panel_hex.rds) produit par
#'     src/covariate_build.R -> lu tel quel, aucune agrégation refaite.
#'   - sinon (données simulées ou cache absent) : construit à la volée depuis les
#'     covariables déclarées dans `data$covariates`.
#'
#' @return Liste : $panel (data.frame long), $hex (grille sf),
#'   $groups (liste nommée covariable -> groupes de modèles concernés).
prepare.modeling.data <- function(cfg, data, hex) {
  cache <- if (!isTRUE(cfg$simulate$use_simulated)) load.covariate.panel(cfg) else NULL

  if (!is.null(cache)) {
    message(">> Covariables : PANEL EN CACHE (", nrow(cache$data), " lignes)")
    groups <- cache$meta$covariate_groups
    return(list(panel = cache$data, hex = cache$hex, groups = groups))
  }

  # Construction à la volée (cas simulé / sans cache).
  years <- sort(unique(st_drop_geometry(data$target)[[data$target_cfg$year_field]]))
  built <- build.hex.panel(hex, data$target, data$target_cfg, data$covariates, years,
                           min_coverage = cfg$hex_grid$min_coverage %||% 0.25)
  # Sans cache, chaque covariable sert aux deux groupes.
  groups <- lapply(names(data$covariates), function(x) c("group1", "group2"))
  names(groups) <- names(data$covariates)
  list(panel = built$data, hex = built$hex, groups = groups)
}

#' Noms des covariables concernées par un groupe de modèles donné
covariates.for.group <- function(groups, which_group) {
  if (length(groups) == 0) return(character(0))
  names(groups)[sapply(groups, function(g) which_group %in% g)]
}


# =============================================================================
# SOUS-PIPELINE GROUPE 1 — prédiction t+1
# =============================================================================
run.group1 <- function(panel, hex, cov_names, cfg) {
  message("\n===== GROUPE 1 : prédiction spatio-temporelle t+1 =====")
  message("  Covariables : ", if (length(cov_names)) paste(cov_names, collapse = ", ") else "(aucune)")

  n_states <- length(levels(panel$target))   # Aucun + niveaux

  # 1. Paires de transition (t -> t+1) avec pression de voisinage.
  pairs <- build.transition.pairs(panel, hex, cov_names)

  # 2. Découpage TEMPOREL (jamais aléatoire).
  splits <- temporal.split(pairs, cfg$models$group1$temporal_split)
  message("  Découpage temporel — train: ", nrow(splits$train),
          " | valid: ", nrow(splits$valid), " | test: ", nrow(splits$test), " paires")

  # 3. Ajustement des modèles demandés.
  run <- cfg$models$group1$run %||% c("markov", "rf", "convlstm")
  results <- list()
  if ("markov"   %in% run) results$markov   <- fit.markov(splits, n_states)$result
  if ("rf"       %in% run) results$rf       <- fit.rf(splits, cov_names, n_states)$result
  if ("convlstm" %in% run) results$convlstm <- fit.convlstm(splits, n_states, hex, cov_names, cfg)$result

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


# =============================================================================
# SOUS-PIPELINE GROUPE 2 — analyse à un instant t
# =============================================================================
run.group2 <- function(panel, hex, cov_names, cfg) {
  message("\n===== GROUPE 2 : analyse économétrique à un instant t =====")
  message("  Covariables : ", if (length(cov_names)) paste(cov_names, collapse = ", ") else "(aucune)")
  year_t <- cfg$models$group2$year_t

  df <- panel[panel$year == year_t, ]
  message("  Année analysée : ", year_t, " (", nrow(df), " hexagones)")
  run <- cfg$models$group2$run %||%
    c("ols","slx","sar","sem","durbin","durbin_error","gam")
  fit.group2(df, hex, cov_names, models = run)
}


# =============================================================================
# EXÉCUTION
# =============================================================================
main <- function(config_path = "config.yaml") {
  cfg <- read_yaml(config_path)

  # 1. Données (réelles ou simulées).
  data <- load.project.data(cfg)

  # 2. Grille hexagonale commune.
  hex <- create.hex.grid(data$zone,
                         cellsize   = cfg$hex_grid$cellsize,
                         target_crs = cfg$project$target_crs %||% 32198)

  # 3. Panel de modélisation (cache réel ou construction à la volée).
  md <- prepare.modeling.data(cfg, data, hex)
  cov1 <- covariates.for.group(md$groups, "group1")   # prédiction (toutes)
  cov2 <- covariates.for.group(md$groups, "group2")   # instant t (sous-ensemble)

  # 4. Sous-pipelines (comparés séparément).
  res1 <- run.group1(md$panel, md$hex, cov1, cfg)
  res2 <- run.group2(md$panel, md$hex, cov2, cfg)

  # 4. Sorties standardisées (par groupe).
  outdir <- cfg$output$dir %||% "data/output/models"
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res1, file.path(outdir, "comparaison_groupe1.csv"), row.names = FALSE)
  write.csv(res2, file.path(outdir, "comparaison_groupe2.csv"), row.names = FALSE)

  message("\n===== RÉSULTATS GROUPE 1 (prédiction t+1) =====")
  print(res1, row.names = FALSE)
  message("\n===== RÉSULTATS GROUPE 2 (instant t) =====")
  print(res2, row.names = FALSE)
  message("\nSorties écrites dans : ", outdir)

  invisible(list(group1 = res1, group2 = res2))
}

# Exécution auto si lancé via Rscript (pas si source() ailleurs).
if (sys.nframe() == 0) {
  main("config.yaml")
}

