# =============================================================================
# src/modelisation.R
# -----------------------------------------------------------------------------
# ORCHESTRATION DE LA MODÉLISATION — deux sous-pipelines comparés séparément :
#
#   GROUPE 2 — ANALYSE à un instant t (inférence) : MCO, SLX, SAR, SEM, Durbin
#              spatial, Durbin erreur, GAM. Métriques : AIC/BIC, pseudo-R²,
#              test de Moran sur les résidus, coefficients spatiaux.
#
#   GROUPE 1 — PRÉDICTION t+1 : Markov spatial, Random Forest, ConvLSTM.
#              Évaluation par validation croisée spatiale + test temporel
#              (kappa pondéré, matrices de confusion, rappel par classe).
#
# Point d'entrée : run.modelisation(). Appelé par main.R dans sa section
# « Modélisation ». Fonctionne sur les VRAIES données (panel en cache produit par
# src/covariate_build.R) ou, à défaut, sur des données SIMULÉES (démo).
#
# Dépendances : sf, yaml + modules hex_grid, features, models_group1/2,
#               model_evaluation (qui charge visualisation + convlstm).
# =============================================================================

library(sf)
library(yaml)

source("src/hex_grid.R")          # grille + agrégation spatiale
source("src/features.R")          # assemblage cible + covariables, cache
source("src/models_group2.R")     # modèles économétriques (instant t)
source("src/model_evaluation.R")  # évaluation groupe 1 (charge models_group1,
                                  # visualisation, model_convlstm)

`%||%` <- function(a, b) if (is.null(a)) b else a


# CHARGEMENT DES DONNÉES (réelles OU simulées) ----------------------------

#' Charger les données du projet selon la config
#'
#' @return Liste : $zone, $target (sf), $covariates, $target_cfg.
load.project.data <- function(cfg) {

  if (isTRUE(cfg$simulate$use_simulated)) {
    # --- Mode DÉMO : données factices ---
    message(">> Mode DONNÉES SIMULÉES (simulate$use_simulated = TRUE)")
    source("src/simulate_data.R")
    sim <- simulate.tbe.dataset(
      seed      = cfg$simulate$seed      %||% 42,
      n_years   = cfg$simulate$n_years   %||% 12,
      extent_km = cfg$simulate$extent_km %||% 200)
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
  # Harmonisation des libellés d'intensité (accents/casse) -> ordinal_levels.
  vf <- cfg$target$value_field
  target[[vf]] <- normalize.intensity(target[[vf]], cfg$target$ordinal_levels)
  target <- st_transform(target, target_crs)

  # Zone d'étude (découpage administratif) ; découpe la cible pour alléger.
  zone <- load.study.zone(cfg$study_zone, target_crs)
  if (is.null(zone)) {
    message("  (study_zone non défini -> zone = emprise de la cible)")
    zone <- st_as_sfc(st_bbox(target)) |> st_sf(geometry = _)
  } else {
    target <- target[zone, ]
    message("  Cible découpée à la zone : ", nrow(target), " entités TBE conservées")
  }
  list(zone = zone, target = target, covariates = cfg$covariates,
       target_cfg = cfg$target)
}


# PRÉPARATION DU PANEL DE MODÉLISATION (hexagone × année) ----------------

#' Obtenir le panel (cible + covariables) : cache réel ou construction à la volée
#'
#' @return Liste : $panel (data.frame long), $hex (grille sf),
#'   $groups (covariable -> groupes de modèles concernés).
prepare.modeling.data <- function(cfg, data, hex) {
  cache <- if (!isTRUE(cfg$simulate$use_simulated)) load.covariate.panel(cfg) else NULL

  if (!is.null(cache)) {
    message(">> Covariables : PANEL EN CACHE (", nrow(cache$data), " lignes)")
    return(list(panel = cache$data, hex = cache$hex, groups = cache$meta$covariate_groups))
  }

  # Sans cache (simulé) : agrégation à la volée ; chaque covariable sert aux 2 groupes.
  years <- sort(unique(st_drop_geometry(data$target)[[data$target_cfg$year_field]]))
  built <- build.hex.panel(hex, data$target, data$target_cfg, data$covariates, years,
                           min_coverage = cfg$hex_grid$min_coverage %||% 0.25)
  groups <- setNames(lapply(names(data$covariates),
                            function(x) c("group1", "group2")), names(data$covariates))
  list(panel = built$data, hex = built$hex, groups = groups)
}

#' Noms des covariables concernées par un groupe de modèles donné
covariates.for.group <- function(groups, which_group) {
  if (length(groups) == 0) return(character(0))
  names(groups)[sapply(groups, function(g) which_group %in% g)]
}


# SOUS-PIPELINE GROUPE 1 — prédiction t+1 (comparaison rapide, un seul split) ----
# NB : l'évaluation APPROFONDIE (CV spatiale + test temporel + graphiques) passe
# par run.model.evaluation() (src/model_evaluation.R). run.group1() ci-dessous
# n'est utilisé qu'en mode simulé (pas de cache), comme comparaison rapide.
run.group1 <- function(panel, hex, cov_names, cfg) {
  message("\n===== GROUPE 1 : prédiction spatio-temporelle t+1 (comparaison rapide) =====")
  n_states <- length(levels(panel$target))
  pairs <- build.transition.pairs(panel, hex, cov_names)
  splits <- temporal.split(pairs, cfg$models$group1$temporal_split)
  message("  Découpage temporel — train: ", nrow(splits$train),
          " | valid: ", nrow(splits$valid), " | test: ", nrow(splits$test), " paires")
  run <- cfg$models$group1$run %||% c("markov", "rf", "convlstm")
  results <- list()
  if ("markov"   %in% run) results$markov   <- fit.markov(splits, n_states)$result
  if ("rf"       %in% run) results$rf       <- fit.rf(splits, cov_names, n_states)$result
  if ("convlstm" %in% run) results$convlstm <- fit.convlstm(splits, n_states, hex, cov_names, cfg)$result
  out <- do.call(rbind, results); rownames(out) <- NULL
  out
}


# SOUS-PIPELINE GROUPE 2 — analyse économétrique à un instant t ---------
run.group2 <- function(panel, hex, cov_names, cfg) {
  message("\n===== GROUPE 2 : analyse économétrique à un instant t =====")
  message("  Covariables : ", if (length(cov_names)) paste(cov_names, collapse = ", ") else "(aucune)")
  year_t <- cfg$models$group2$year_t
  df <- panel[panel$year == year_t, ]
  message("  Année analysée : ", year_t, " (", nrow(df), " hexagones)")
  run <- cfg$models$group2$run %||% c("ols","slx","sar","sem","durbin","durbin_error","gam")
  fit.group2(df, hex, cov_names, models = run)
}


# ORCHESTRATION COMPLÈTE ----------------------------------------------------

#' Lancer toute la modélisation (groupe 2 + groupe 1) et écrire les sorties
#'
#' Sorties (sous cfg$output$root) :
#'   economie/    -> comparaison des modèles économétriques (groupe 2)
#'   prediction/  -> évaluation des modèles prédictifs (groupe 1) + graphiques
#'   cartes/      -> cartes territoriales observé/prédit/erreur
#'
#' @param config_path Chemin du YAML.
run.modelisation <- function(config_path = "config.yaml") {
  cfg  <- read_yaml(config_path)
  root <- cfg$output$root %||% "outputs"

  # --- Données + grille + panel ---
  data <- load.project.data(cfg)
  hex  <- create.hex.grid(data$zone, cfg$hex_grid$cellsize,
                          cfg$project$target_crs %||% 32198)
  md   <- prepare.modeling.data(cfg, data, hex)

  # --- GROUPE 2 : économétrie (instant t) ---
  cov2 <- covariates.for.group(md$groups, "group2")
  res2 <- run.group2(md$panel, md$hex, cov2, cfg)
  dir.create(file.path(root, "economie"), recursive = TRUE, showWarnings = FALSE)
  write.csv(res2, file.path(root, "economie", "comparaison_groupe2.csv"), row.names = FALSE)
  message("\n== GROUPE 2 (instant t) =="); print(res2, row.names = FALSE)

  # --- GROUPE 1 : prédiction (évaluation approfondie si cache, sinon rapide) ---
  cache <- if (!isTRUE(cfg$simulate$use_simulated)) load.covariate.panel(cfg) else NULL
  if (!is.null(cache)) {
    run.model.evaluation(config_path)     # CV + test temporel + graphiques -> prediction/
    map.territory.results(config_path)    # cartes observé/prédit/erreur -> cartes/
  } else {
    cov1 <- covariates.for.group(md$groups, "group1")
    res1 <- run.group1(md$panel, md$hex, cov1, cfg)
    dir.create(file.path(root, "prediction"), recursive = TRUE, showWarnings = FALSE)
    write.csv(res1, file.path(root, "prediction", "comparaison_groupe1.csv"), row.names = FALSE)
    message("\n== GROUPE 1 (prédiction, rapide) =="); print(res1, row.names = FALSE)
  }

  message("\nModélisation terminée. Sorties sous : ", root, "/")
  invisible(TRUE)
}

# Exécution directe (Rscript src/modelisation.R) : lance toute la modélisation.
if (sys.nframe() == 0) run.modelisation("config.yaml")
