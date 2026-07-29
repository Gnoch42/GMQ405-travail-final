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
  target <- st_read(cfg$target$source, quiet = TRUE)
  # Normalisation des libellés d'intensité (accents/casse incohérents dans la
  # source : "Leger"/"Léger", "Modere"/"Modéré"). On harmonise pour l'ordinal.
  vf <- cfg$target$value_field
  target[[vf]] <- normalize.intensity(target[[vf]], cfg$target$ordinal_levels)

  # La zone d'étude par défaut = emprise de la cible (les covariables seront
  # découpées dessus lors de l'agrégation).
  zone <- st_as_sfc(st_bbox(target)) |> st_sf(geometry = _)

  list(zone = zone, target = target, covariates = cfg$covariates,
       target_cfg = cfg$target)
}


#' Harmoniser les libellés d'intensité (gère accents et casse variables)
normalize.intensity <- function(x, levels) {
  x <- as.character(x)
  # Table de correspondance : variantes sans accent -> libellé canonique.
  key <- gsub("é", "e", tolower(trimws(x)))
  canon <- setNames(levels, gsub("é", "e", tolower(levels)))
  out <- canon[key]
  out[is.na(out)] <- x[is.na(out)]   # on laisse tel quel si non reconnu
  unname(out)
}

# Opérateur "valeur par défaut si NULL" (pratique pour lire le YAML).
`%||%` <- function(a, b) if (is.null(a)) b else a


# =============================================================================
# SOUS-PIPELINE GROUPE 1 — prédiction t+1
# =============================================================================
run.group1 <- function(hex, data, cfg) {
  message("\n===== GROUPE 1 : prédiction spatio-temporelle t+1 =====")
  cov_names <- names(data$covariates)

  # 1. Panneau long hexagone×année (cible + covariables).
  years <- sort(unique(st_drop_geometry(data$target)[[data$target_cfg$year_field]]))
  panel <- build.hex.panel(hex, data$target, data$target_cfg,
                           data$covariates, years,
                           min_coverage = cfg$hex_grid$min_coverage %||% 0.25)

  n_states <- length(levels(panel$data$target))   # Aucun + niveaux

  # 2. Paires de transition (t -> t+1) avec pression de voisinage.
  pairs <- build.transition.pairs(panel$data, panel$hex, cov_names)

  # 3. Découpage TEMPOREL (jamais aléatoire).
  splits <- temporal.split(pairs, cfg$models$group1$temporal_split)
  message("  Découpage temporel — train: ", nrow(splits$train),
          " | valid: ", nrow(splits$valid), " | test: ", nrow(splits$test), " paires")

  # 4. Ajustement des modèles demandés.
  run <- cfg$models$group1$run %||% c("markov", "rf", "convlstm")
  results <- list()
  if ("markov"   %in% run) results$markov   <- fit.markov(splits, n_states)$result
  if ("rf"       %in% run) results$rf       <- fit.rf(splits, cov_names, n_states)$result
  if ("convlstm" %in% run) results$convlstm <- fit.convlstm(splits, n_states)$result

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


# =============================================================================
# SOUS-PIPELINE GROUPE 2 — analyse à un instant t
# =============================================================================
run.group2 <- function(hex, data, cfg) {
  message("\n===== GROUPE 2 : analyse économétrique à un instant t =====")
  cov_names <- names(data$covariates)
  year_t <- cfg$models$group2$year_t

  # Panneau restreint à l'année d'analyse.
  panel <- build.hex.panel(hex, data$target, data$target_cfg,
                           data$covariates, year_t,
                           min_coverage = cfg$hex_grid$min_coverage %||% 0.25)
  df <- panel$data[panel$data$year == year_t, ]

  message("  Année analysée : ", year_t, " (", nrow(df), " hexagones)")
  run <- cfg$models$group2$run %||%
    c("ols","slx","sar","sem","durbin","durbin_error","gam")
  fit.group2(df, panel$hex, cov_names, models = run)
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

  # 3. Sous-pipelines (comparés séparément).
  res1 <- run.group1(hex, data, cfg)
  res2 <- run.group2(hex, data, cfg)

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
