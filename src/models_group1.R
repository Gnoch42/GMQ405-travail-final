# =============================================================================
# src/models_group1.R
# -----------------------------------------------------------------------------
# GROUPE 1 — Modèles de PRÉDICTION spatio-temporelle t+1.
#
# But de ce groupe : PRÉDIRE où sera la TBE l'année suivante. On apprend la
# dynamique passée (année t -> année t+1) et on valide sur des années futures
# jamais vues à l'entraînement.
#
# Les 3 modèles :
#   MARKOV (CA-Markov) : matrice de probabilités de transition entre classes
#                        d'intensité (t -> t+1), modulée par la pression de
#                        voisinage. Simple et interprétable.
#   RF / XGBOOST       : apprentissage tabulaire ; chaque ligne = hexagone×année,
#                        cible = intensité en t+1. Gère les interactions non
#                        linéaires. (Nécessite le package ranger ou xgboost.)
#   ConvLSTM           : réseau de neurones sur séquences de cartes rasterisées.
#                        (Nécessite keras3/torch + rasterisation ; fourni en
#                        SQUELETTE, voir la garde plus bas.)
#
# ---------------------------------------------------------------------------
# POINT MÉTHODOLOGIQUE CENTRAL : le DÉCOUPAGE TEMPOREL (jamais aléatoire)
# ---------------------------------------------------------------------------
# On prédit le futur : l'entraînement, la validation et le test doivent donc
# être séparés DANS LE TEMPS (blocs d'années consécutives). Un découpage
# aléatoire mélangerait passé et futur et surestimerait gravement la performance.
# Voir temporal.split() ci-dessous.
#
# Dépendances : sf, spdep ; ranger/xgboost/keras3 optionnels (gardes prévues).
# =============================================================================

library(sf)


# -----------------------------------------------------------------------------
# DÉCOUPAGE TEMPOREL
# -----------------------------------------------------------------------------

#' Séparer les paires (t -> t+1) en entraînement / validation / test par années
#'
#' @param pairs  data.frame de paires (voir build.transition.pairs()), colonne
#'               `year_t` = année de départ de la transition.
#' @param split  Liste : train_years=c(min,max), valid_years=c(min,max),
#'               test_years=c(min,max) — bornes INCLUSIVES sur year_t.
#' @return La liste des trois sous-ensembles (train, valid, test).
temporal.split <- function(pairs, split) {
  in_range <- function(y, rng) y >= rng[1] & y <= rng[2]
  list(
    train = pairs[in_range(pairs$year_t, split$train_years), ],
    valid = pairs[in_range(pairs$year_t, split$valid_years), ],
    test  = pairs[in_range(pairs$year_t, split$test_years), ]
  )
}


# -----------------------------------------------------------------------------
# CONSTRUCTION DES PAIRES DE TRANSITION (t -> t+1)
# -----------------------------------------------------------------------------

#' Construire les paires hexagone (année t, année t+1) avec voisinage
#'
#' Transforme le panneau long (hexagone×année) en un tableau où chaque ligne
#' contient l'état en t, l'état futur en t+1 (la cible à prédire), et un indice
#' de PRESSION DE VOISINAGE en t (part de voisins infestés) — variable clé pour
#' la propagation spatiale.
#'
#' @param panel data.frame issu de build.hex.panel() (hex_id, year, target, cov).
#' @param hex   Grille (géométrie) pour calculer le voisinage.
#' @param cov_names Covariables à reporter comme variables explicatives.
#' @return data.frame : hex_id, year_t, state_t, state_t1 (cible), neigh_pressure,
#'   + covariables.
build.transition.pairs <- function(panel, hex, cov_names) {
  # Voisinage (liste des voisins de chaque hexagone).
  nb <- spdep::poly2nb(hex, queen = TRUE)
  names(nb) <- as.character(hex$hex_id)

  panel$state <- as.integer(panel$target) - 1L    # 0 = Aucun, 1..k
  years <- sort(unique(panel$year))

  rows <- list()
  for (i in seq_len(length(years) - 1)) {
    yt <- years[i]; yt1 <- years[i + 1]
    a <- panel[panel$year == yt,  ]
    b <- panel[panel$year == yt1, ]
    # Aligner b sur a par hex_id (même ordre).
    b <- b[match(a$hex_id, b$hex_id), ]

    # Pression de voisinage en t : proportion de voisins infestés (state > 0).
    infested <- setNames(a$state > 0, as.character(a$hex_id))
    pressure <- sapply(seq_len(nrow(a)), function(j) {
      hid <- as.character(a$hex_id[j])
      vois <- nb[[hid]]
      if (length(vois) == 0 || vois[1] == 0) return(0)
      mean(infested[as.character(hex$hex_id[vois])], na.rm = TRUE)
    })

    df <- data.frame(
      hex_id   = a$hex_id,
      year_t   = yt,
      state_t  = a$state,
      state_t1 = b$state,               # CIBLE : intensité l'année suivante
      neigh_pressure = round(pressure, 3)
    )
    for (cn in cov_names) df[[cn]] <- a[[cn]]
    rows[[as.character(yt)]] <- df
  }
  out <- do.call(rbind, rows)
  out[!is.na(out$state_t1), ]           # on retire les transitions incomplètes
}


# -----------------------------------------------------------------------------
# MODÈLE 1 — MARKOV SPATIAL (CA-Markov)
# -----------------------------------------------------------------------------

#' Ajuster une matrice de transition de Markov et prédire t+1
#'
#' Principe : on compte, sur l'entraînement, combien de fois un hexagone passe
#' de l'état i (en t) à l'état j (en t+1). Normalisé par ligne, cela donne la
#' PROBABILITÉ de transition i -> j. La prédiction retient l'état le plus
#' probable (MAP). Version simple ; la pression de voisinage pourrait moduler
#' ces probabilités (extension via régression ordinale, cf. méthodologie).
#'
#' @return Liste : $matrix (matrice de transition), $result (ligne standardisée).
fit.markov <- function(splits, n_states) {
  train <- splits$train
  states <- 0:(n_states - 1)

  # Matrice de comptage des transitions observées.
  trans <- table(factor(train$state_t,  levels = states),
                 factor(train$state_t1, levels = states))
  # Normalisation par ligne -> probabilités (gère les lignes vides).
  probs <- prop.table(trans + 1e-9, margin = 1)   # +epsilon évite division par 0

  # Prédiction sur le test : état le plus probable depuis l'état de départ.
  test <- splits$test
  pred <- states[apply(probs[as.character(test$state_t), , drop = FALSE], 1, which.max)]

  res <- standardize.group1.result("Markov spatial", test$state_t1, pred, n_states)
  list(matrix = probs, result = res)
}


# -----------------------------------------------------------------------------
# MODÈLE 2 — RANDOM FOREST / XGBOOST
# -----------------------------------------------------------------------------

#' Ajuster un Random Forest (ou XGBoost) pour prédire l'état t+1
#'
#' Variables explicatives : état en t, pression de voisinage, covariables.
#' Utilise `ranger` si disponible, sinon renvoie un résultat "package_manquant"
#' (le pipeline continue sans planter — comportement voulu pour un squelette).
fit.rf <- function(splits, cov_names, n_states, engine = "ranger") {
  if (!requireNamespace(engine, quietly = TRUE)) {
    message("  [RF] package '", engine, "' absent -> modèle ignoré (installez-le pour l'activer).")
    return(list(result = standardize.group1.result("Random Forest", NULL, NULL,
                                                    n_states, statut = "package_manquant")))
  }

  predictors <- c("state_t", "neigh_pressure", cov_names)
  train <- splits$train; test <- splits$test
  train$state_t1 <- factor(train$state_t1, levels = 0:(n_states - 1))

  fml <- as.formula(paste("state_t1 ~", paste(predictors, collapse = " + ")))
  m <- ranger::ranger(fml, data = train[complete.cases(train[, predictors]), ],
                      num.trees = 300, probability = FALSE)
  pred <- as.integer(as.character(predict(m, test)$predictions))

  res <- standardize.group1.result("Random Forest", test$state_t1, pred, n_states)
  list(model = m, result = res)
}


# -----------------------------------------------------------------------------
# MODÈLE 3 — ConvLSTM (SQUELETTE)
# -----------------------------------------------------------------------------

#' Squelette ConvLSTM
#'
#' Le ConvLSTM travaille sur des SÉQUENCES DE CARTES (tenseur temps × lignes ×
#' colonnes × canaux) et nécessite : (1) une rasterisation régulière des
#' hexagones, (2) keras3/torch. Ces dépendances lourdes ne sont pas installées
#' par défaut. On fournit ici la STRUCTURE et un garde-fou : la fonction renvoie
#' un résultat standardisé "non_implemente" pour que la comparaison de groupe 1
#' reste cohérente sans bloquer le pipeline.
#'
#' Étapes à implémenter (voir méthodologie) :
#'   1. Rasteriser la cible par année en cartes de mêmes dimensions.
#'   2. Empiler en tenseur (t, h, w, canaux) + canaux covariables.
#'   3. Découpage temporel strict (déjà géré par temporal.split en amont).
#'   4. Définir l'architecture (couches ConvLSTM2D) et entraîner.
#'   5. Prédire la carte t+1, reconvertir en hexagones, évaluer (kappa pondéré).
fit.convlstm <- function(splits, n_states) {
  has_keras <- requireNamespace("keras3", quietly = TRUE) ||
               requireNamespace("torch",  quietly = TRUE)
  if (!has_keras) {
    message("  [ConvLSTM] keras3/torch absent -> squelette non exécuté (voir doc).")
    return(list(result = standardize.group1.result("ConvLSTM", NULL, NULL,
                                                    n_states, statut = "non_implemente")))
  }
  # (Implémentation deep learning à compléter ici.)
  list(result = standardize.group1.result("ConvLSTM", NULL, NULL,
                                           n_states, statut = "a_completer"))
}


# -----------------------------------------------------------------------------
# STANDARDISATION DE LA SORTIE (structure commune au GROUPE 1)
# -----------------------------------------------------------------------------
# Métriques adaptées à une cible ORDINALE, comparables entre les 3 modèles :
#   - accuracy       : proportion d'hexagones dont la classe est exactement bonne.
#   - kappa_pondere  : accord observé/attendu au-delà du hasard, PÉNALISANT
#                      davantage les grosses erreurs (Aucun prédit "Grave") que
#                      les petites (Léger prédit "Modéré"). Métrique de référence
#                      pour l'ordinal. 1 = parfait, 0 = hasard.
#   - rmse_rang      : erreur quadratique moyenne sur les RANGS d'intensité.

#' Kappa pondéré quadratique (sans dépendance externe)
#'
#' Calcule l'accord corrigé du hasard entre observé et prédit, avec des poids
#' quadratiques : une erreur de 2 classes compte 4× une erreur de 1 classe.
weighted.kappa <- function(obs, pred, n_states) {
  lv <- 0:(n_states - 1)
  O <- table(factor(obs, levels = lv), factor(pred, levels = lv))
  N <- sum(O)
  if (N == 0) return(NA_real_)
  O <- O / N
  # Matrice attendue sous indépendance (produit des marges).
  E <- outer(rowSums(O), colSums(O))
  # Poids quadratiques normalisés.
  W <- outer(lv, lv, function(i, j) (i - j)^2) / (n_states - 1)^2
  1 - sum(W * O) / sum(W * E)
}

standardize.group1.result <- function(name, obs, pred, n_states, statut = "ok") {
  if (is.null(obs) || is.null(pred) || length(obs) == 0) {
    return(data.frame(modele = name, statut = statut,
                      accuracy = NA_real_, kappa_pondere = NA_real_,
                      rmse_rang = NA_real_, n_test = 0L, stringsAsFactors = FALSE))
  }
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]; pred <- pred[ok]
  data.frame(
    modele        = name,
    statut        = statut,
    accuracy      = round(mean(obs == pred), 3),
    kappa_pondere = round(weighted.kappa(obs, pred, n_states), 3),
    rmse_rang     = round(sqrt(mean((obs - pred)^2)), 3),
    n_test        = length(obs),
    stringsAsFactors = FALSE
  )
}
