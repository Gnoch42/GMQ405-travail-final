# =============================================================================
# src/model_evaluation.R
# -----------------------------------------------------------------------------
# ÉVALUATION comparable des trois modèles de prédiction (groupe 1) sur la cible
# ordinale à 4 états : Absence < Légère < Modérée < Sévère.
#
# Fournit des fonctions RÉUTILISABLES (appelables avec observations + prédictions)
# — à relancer chaque fois que les modèles sont ajustés — pour :
#   1. Matrice de confusion (comptes bruts + normalisée par ligne, en %).
#   2. Kappa pondéré QUADRATIQUE (métrique principale, car classes ordinales).
#   3. Précision et rappel PAR CLASSE (détecte les faiblesses sur classes rares).
#   4. Kappa simulation (décomposition quantité / localisation, Pontius 2000) —
#      pour le modèle Markov (CA-Markov).
#
# Et une VALIDATION CROISÉE SPATIO-TEMPORELLE identique pour les trois modèles :
#   - CV par BLOCS SPATIAUX contigus (évite la fuite d'information due à
#     l'autocorrélation entre hexagones voisins) -> distribution de performance.
#   - TEST TEMPOREL final sur la phase de déclin du cycle (années réservées) ->
#     mesure la vraie capacité prédictive de la dynamique, pas l'interpolation.
#
# Dépendances : sf, ggplot2 ; + hex_grid.R, features.R, models_group1.R.
# =============================================================================

library(sf)
library(ggplot2)

source("src/hex_grid.R")
source("src/features.R")
source("src/models_group1.R")   # interfaces markov.predict / rf.predict / convlstm.predict
try(source("src/model_convlstm.R"), silent = TRUE)   # ConvLSTM (torch), si disponible

`%||%` <- function(a, b) if (is.null(a)) b else a

# Libellés des 4 états ordonnés (index 0..3).
ETATS_TBE <- c("Absence", "Légère", "Modérée", "Sévère")


# =============================================================================
# MÉTRIQUES ÉLÉMENTAIRES (appelables avec obs + pred)
# =============================================================================

#' Matrice de confusion : comptes bruts + version normalisée par ligne (%)
#'
#' Lignes = état OBSERVÉ, colonnes = état PRÉDIT. La normalisation par ligne
#' donne, pour chaque état réel, la répartition (%) des prédictions : la diagonale
#' est le rappel de la classe.
#'
#' @param obs,pred Vecteurs d'états entiers (0..n_states-1).
#' @param n_states Nombre d'états.
#' @param labels   Libellés des états (défaut ETATS_TBE).
#' @return Liste : $counts (matrice), $row_pct (matrice %).
confusion.matrix <- function(obs, pred, n_states, labels = ETATS_TBE) {
  lv <- 0:(n_states - 1)
  ok <- !is.na(obs) & !is.na(pred)
  cm <- table(Observe = factor(obs[ok], levels = lv, labels = labels[1:n_states]),
              Predit  = factor(pred[ok], levels = lv, labels = labels[1:n_states]))
  row_pct <- round(100 * prop.table(cm, margin = 1), 1)
  row_pct[is.nan(row_pct)] <- 0
  list(counts = cm, row_pct = row_pct)
}


#' Kappa de Cohen PONDÉRÉ (pondération quadratique) — métrique principale
#'
#' Corrige l'accord du hasard et pénalise DAVANTAGE les grosses erreurs
#' (prédire "Absence" là où c'est "Sévère") que les petites (Modérée vs Sévère).
#' Poids quadratiques : une erreur de d classes coûte d². 1 = parfait, 0 = hasard,
#' négatif = pire que le hasard.
kappa.weighted.quad <- function(obs, pred, n_states) {
  lv <- 0:(n_states - 1)
  ok <- !is.na(obs) & !is.na(pred)
  if (sum(ok) == 0) return(NA_real_)
  O <- table(factor(obs[ok], levels = lv), factor(pred[ok], levels = lv))
  N <- sum(O); O <- O / N
  E <- outer(rowSums(O), colSums(O))               # accord attendu (marges)
  W <- outer(lv, lv, function(i, j) (i - j)^2) / (n_states - 1)^2   # poids quadratiques
  denom <- sum(W * E)
  if (denom == 0) return(NA_real_)
  1 - sum(W * O) / denom
}


#' Précision, rappel, F1 et effectif PAR CLASSE
#'
#' - Rappel (sensibilité) = part des cas réels de la classe correctement prédits.
#' - Précision = part des prédictions de la classe qui sont correctes.
#' Utile pour repérer une sous-performance sur les classes rares (ex. "Sévère").
class.metrics <- function(obs, pred, n_states, labels = ETATS_TBE) {
  lv <- 0:(n_states - 1)
  ok <- !is.na(obs) & !is.na(pred); obs <- obs[ok]; pred <- pred[ok]
  rows <- lapply(lv, function(k) {
    tp <- sum(obs == k & pred == k)
    fn <- sum(obs == k & pred != k)
    fp <- sum(obs != k & pred == k)
    prec <- if (tp + fp == 0) NA_real_ else tp / (tp + fp)
    rec  <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
    f1   <- if (is.na(prec) || is.na(rec) || prec + rec == 0) NA_real_ else 2 * prec * rec / (prec + rec)
    data.frame(classe = labels[k + 1], effectif = tp + fn,
               precision = round(prec, 3), rappel = round(rec, 3), F1 = round(f1, 3))
  })
  do.call(rbind, rows)
}


#' Kappa simulation — décomposition QUANTITÉ / LOCALISATION (Pontius 2000)
#'
#' Standard pour valider les modèles CA-Markov : distingue deux sources d'erreur.
#'   - Kquantity : le modèle prédit-il la BONNE PROPORTION de chaque état ?
#'                 (bon nombre d'hexagones "Sévère", indépendamment de l'endroit)
#'   - Klocation : place-t-il ces états au BON ENDROIT, une fois la quantité fixée ?
#' Notations (proportions de la matrice de confusion) :
#'   P(p) = accord observé (diagonale) ; P(k) = accord attendu au hasard (produit
#'   des marges) ; P(m) = accord maximal possible avec les marges = Σ min(marges).
#'   Klocation = (P(p) − P(k)) / (P(m) − P(k)) ; Kquantity = (P(m) − P(k))/(1 − P(k)).
kappa.simulation <- function(obs, pred, n_states) {
  lv <- 0:(n_states - 1)
  ok <- !is.na(obs) & !is.na(pred)
  if (sum(ok) == 0) return(list(kappa = NA, k_location = NA, k_quantity = NA))
  O <- table(factor(obs[ok], levels = lv), factor(pred[ok], levels = lv))
  P <- O / sum(O)
  Pp <- sum(diag(P))                       # accord observé
  r <- rowSums(P); s <- colSums(P)
  Pk <- sum(r * s)                         # accord attendu (hasard)
  Pm <- sum(pmin(r, s))                    # accord max donné les marges
  klocation <- if (Pm - Pk == 0) NA_real_ else (Pp - Pk) / (Pm - Pk)
  kquantity <- if (1 - Pk == 0)  NA_real_ else (Pm - Pk) / (1 - Pk)
  kstandard <- if (1 - Pk == 0)  NA_real_ else (Pp - Pk) / (1 - Pk)
  list(kappa = round(kstandard, 3),
       k_location = round(klocation, 3), k_quantity = round(kquantity, 3))
}


#' Évaluation complète d'un jeu de prédictions (toutes les métriques ci-dessus)
#'
#' @param model_name Nom du modèle (pour l'étiquetage).
#' @param pontius    TRUE pour ajouter le Kappa simulation (modèle Markov).
#' @return Liste : $summary (1 ligne), $confusion, $class, [$pontius].
evaluate.predictions <- function(obs, pred, n_states, model_name, pontius = FALSE,
                                 labels = ETATS_TBE) {
  if (is.null(pred) || length(pred) == 0 || all(is.na(pred))) {
    return(list(summary = data.frame(modele = model_name, statut = "non_disponible",
                                     kappa_pondere = NA, accuracy = NA),
                confusion = NULL, class = NULL))
  }
  cm  <- confusion.matrix(obs, pred, n_states, labels)
  kap <- kappa.weighted.quad(obs, pred, n_states)
  cls <- class.metrics(obs, pred, n_states, labels)
  ok  <- !is.na(obs) & !is.na(pred)
  acc <- mean(obs[ok] == pred[ok])
  out <- list(summary = data.frame(modele = model_name, statut = "ok",
                                   kappa_pondere = round(kap, 3), accuracy = round(acc, 3)),
              confusion = cm, class = cls)
  if (pontius) out$pontius <- kappa.simulation(obs, pred, n_states)
  out
}


# =============================================================================
# DÉCOUPAGE SPATIAL EN BLOCS CONTIGUS
# =============================================================================

#' Partitionner les hexagones en k blocs spatiaux CONTIGUS
#'
#' On regroupe les hexagones par proximité géographique (k-means sur les
#' coordonnées des centroïdes) : chaque bloc est une zone d'un seul tenant. C'est
#' préférable à un tirage aléatoire d'hexagones, qui laisserait des voisins
#' corrélés à la fois dans l'entraînement et le test (fuite d'information due à
#' l'autocorrélation spatiale) et gonflerait artificiellement la performance.
#'
#' @return Vecteur nommé hex_id -> numéro de bloc (1..k).
make.spatial.blocks <- function(hex, k = 4, seed = 1) {
  set.seed(seed)
  xy <- st_coordinates(st_centroid(st_geometry(hex)))
  k <- min(k, nrow(unique(xy)))
  cl <- kmeans(xy, centers = k, nstart = 10)$cluster
  setNames(cl, hex$hex_id)
}


# =============================================================================
# VALIDATION CROISÉE SPATIO-TEMPORELLE
# =============================================================================

#' Registre des modèles à évaluer (interfaces de prédiction communes)
#'
#' Chaque entrée : $name, $predict(train, test) -> vecteur d'états prédits (ou
#' NULL), $pontius (TRUE pour Markov). Le MÊME découpage est appliqué à tous.
build.model.registry <- function(cov_names, n_states, hex = NULL, cfg = NULL) {
  # ConvLSTM : implémentation torch si disponible, sinon renvoie NULL (NA).
  convlstm_fn <- if (exists("convlstm.predict.full") && !is.null(hex)) {
    function(tr, te) convlstm.predict.full(tr, te, n_states, hex, cov_names, cfg)
  } else {
    function(tr, te) NULL
  }
  list(
    list(name = "Markov spatial", pontius = TRUE,
         predict = function(tr, te) markov.predict(tr, te, n_states)),
    list(name = "Random Forest", pontius = FALSE,
         predict = function(tr, te) rf.predict(tr, te, cov_names, n_states)),
    list(name = "ConvLSTM", pontius = FALSE, predict = convlstm_fn)
  )
}


#' Validation croisée spatio-temporelle des trois modèles
#'
#' @param pairs     Paires de transition (build.transition.pairs) avec colonne hex_id.
#' @param hex       Grille (géométrie) pour les blocs spatiaux.
#' @param cov_names Covariables explicatives.
#' @param n_states  Nombre d'états.
#' @param train_years,test_years Bornes temporelles (year_t). test_years = phase
#'   de déclin réservée.
#' @param n_blocks  Nombre de blocs spatiaux pour la CV.
#' @param outdir    Dossier de sortie (tables + images).
#' @return Liste : $summary (tableau récap), $per_block, $confusions, $class, $pontius.
spatiotemporal.evaluation <- function(pairs, hex, cov_names, n_states,
                                      train_years, test_years, n_blocks = 4,
                                      outdir = "data/output/evaluation",
                                      labels = ETATS_TBE, cfg = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  models <- build.model.registry(cov_names, n_states, hex, cfg)

  # Blocs spatiaux + partition temporelle (identiques pour tous les modèles).
  blocks <- make.spatial.blocks(hex, n_blocks)
  pairs$block <- blocks[as.character(pairs$hex_id)]
  in_train <- pairs$year_t >= train_years[1] & pairs$year_t <= train_years[2]
  in_test  <- pairs$year_t >= test_years[1]  & pairs$year_t <= test_years[2]
  train_all <- pairs[in_train, ]
  test_temporal <- pairs[in_test, ]
  ublocks <- sort(unique(train_all$block))
  message("Évaluation : ", nrow(train_all), " paires d'entraînement (",
          length(ublocks), " blocs), ", nrow(test_temporal), " paires de test temporel.")

  summ <- list(); per_block <- list(); confs <- list(); classes <- list(); pontius <- list()

  for (mdl in models) {
    # --- (a) CV par blocs spatiaux (dans la période d'entraînement) ---
    kappas <- sapply(ublocks, function(b) {
      tr <- train_all[train_all$block != b, ]
      te <- train_all[train_all$block == b, ]
      pr <- mdl$predict(tr, te)
      if (is.null(pr)) NA_real_ else kappa.weighted.quad(te$state_t1, pr, n_states)
    })
    per_block[[mdl$name]] <- data.frame(modele = mdl$name, bloc = ublocks,
                                        kappa_pondere = round(kappas, 3))

    # --- (b) Test temporel final (entraîné sur toute la période d'entraînement) ---
    pr_t <- if (nrow(test_temporal) > 0) mdl$predict(train_all, test_temporal) else NULL
    ev <- evaluate.predictions(test_temporal$state_t1, pr_t, n_states, mdl$name,
                               pontius = mdl$pontius, labels = labels)

    # --- Récapitulatif par modèle ---
    k_ok <- kappas[!is.na(kappas)]   # NA propre (et non NaN) si aucun bloc évalué
    summ[[mdl$name]] <- data.frame(
      modele = mdl$name,
      kappa_spatial_moy = if (length(k_ok) == 0) NA_real_ else round(mean(k_ok), 3),
      kappa_spatial_sd  = if (length(k_ok) < 2)  NA_real_ else round(stats::sd(k_ok), 3),
      kappa_temporel    = ev$summary$kappa_pondere,
      accuracy_temporel = ev$summary$accuracy,
      statut = ev$summary$statut)
    confs[[mdl$name]]  <- ev$confusion
    classes[[mdl$name]] <- if (!is.null(ev$class)) cbind(modele = mdl$name, ev$class) else NULL
    if (!is.null(ev$pontius)) pontius[[mdl$name]] <- ev$pontius

    # Sauvegardes par modèle.
    if (!is.null(ev$confusion)) {
      save.confusion(ev$confusion, mdl$name, outdir)
    }
  }

  summary_tab <- do.call(rbind, summ);   rownames(summary_tab) <- NULL
  class_tab   <- do.call(rbind, classes); rownames(class_tab) <- NULL

  # Écritures récapitulatives + graphiques de synthèse.
  write.csv(summary_tab, file.path(outdir, "comparaison_modeles.csv"), row.names = FALSE)
  plot.kappa.comparison(summary_tab, outdir)
  if (!is.null(class_tab)) {
    write.csv(class_tab, file.path(outdir, "precision_rappel_par_classe.csv"), row.names = FALSE)
    plot.class.recall(class_tab, outdir, labels = labels)
  }
  if (length(pontius) > 0) {
    pdf_tab <- do.call(rbind, lapply(names(pontius), function(nm)
      data.frame(modele = nm, as.data.frame(pontius[[nm]]))))
    write.csv(pdf_tab, file.path(outdir, "kappa_simulation_markov.csv"), row.names = FALSE)
  }

  message("\n== Comparaison des modèles =="); print(summary_tab, row.names = FALSE)
  if (!is.null(class_tab)) { message("\n== Précision / rappel par classe (test temporel) =="); print(class_tab, row.names = FALSE) }
  if (length(pontius) > 0) { message("\n== Kappa simulation (Markov) =="); print(pontius) }
  message("\nSorties écrites dans : ", outdir)

  invisible(list(summary = summary_tab, per_block = per_block, confusions = confs,
                 class = class_tab, pontius = pontius))
}


#' Sauvegarder une matrice de confusion en table CSV + image (heatmap)
save.confusion <- function(cm, model_name, outdir) {
  slug <- gsub("[^A-Za-z0-9]+", "_", tolower(model_name))
  write.csv(as.data.frame.matrix(cm$counts),
            file.path(outdir, paste0("confusion_", slug, ".csv")))

  df <- as.data.frame(cm$row_pct)
  names(df) <- c("Observe", "Predit", "pct")
  df$compte <- as.data.frame(cm$counts)$Freq
  p <- ggplot(df, aes(Predit, Observe, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(compte, "\n", pct, "%")), size = 3) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 100)) +
    scale_y_discrete(limits = rev) +
    labs(title = paste("Matrice de confusion —", model_name),
         subtitle = "ligne = observé, colonne = prédit ; % normalisé par ligne (rappel en diagonale)",
         fill = "% ligne") +
    theme_minimal(base_size = 11)
  ggsave(file.path(outdir, paste0("confusion_", slug, ".png")), p,
         width = 6.5, height = 5.5, dpi = 120)
}


#' Graphique de comparaison du KAPPA PONDÉRÉ entre modèles
#'
#' Barres groupées : kappa en validation croisée spatiale (moyenne ± écart-type
#' sur les blocs) et kappa au test temporel (phase de déclin). Métrique principale
#' de comparaison des modèles.
plot.kappa.comparison <- function(summary_tab, outdir) {
  df <- rbind(
    data.frame(modele = summary_tab$modele, type = "CV spatiale",
               kappa = summary_tab$kappa_spatial_moy, sd = summary_tab$kappa_spatial_sd),
    data.frame(modele = summary_tab$modele, type = "Test temporel",
               kappa = summary_tab$kappa_temporel, sd = NA_real_))
  df$modele <- factor(df$modele, levels = summary_tab$modele)
  p <- ggplot(df, aes(modele, kappa, fill = type)) +
    geom_col(position = position_dodge(0.8), width = 0.7, na.rm = TRUE) +
    geom_errorbar(aes(ymin = kappa - sd, ymax = kappa + sd),
                  position = position_dodge(0.8), width = 0.2, na.rm = TRUE) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_text(aes(label = ifelse(is.na(kappa), "n/d", sprintf("%.2f", kappa))),
              position = position_dodge(0.8), vjust = -0.5, size = 3, na.rm = TRUE) +
    scale_fill_manual(values = c("CV spatiale" = "#6baed6", "Test temporel" = "#08519c")) +
    labs(title = "Kappa pondéré par modèle",
         subtitle = "CV spatiale (moyenne ± écart-type sur blocs) vs test temporel (phase de déclin)",
         x = NULL, y = "Kappa pondéré (quadratique)", fill = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  ggsave(file.path(outdir, "kappa_comparaison.png"), p, width = 7.5, height = 5, dpi = 120)
}


#' Graphique du RAPPEL PAR CLASSE (détection des classes rares)
#'
#' Complète les matrices de confusion : montre, pour chaque modèle, la part des
#' cas réels de chaque état correctement détectés. Met en évidence les modèles
#' qui « ratent » les classes rares (Modéré/Grave).
plot.class.recall <- function(class_tab, outdir, labels = ETATS_TBE) {
  df <- class_tab
  df$classe <- factor(df$classe, levels = labels)
  df$modele <- factor(df$modele, levels = unique(df$modele))
  p <- ggplot(df, aes(classe, rappel, fill = modele)) +
    geom_col(position = position_dodge(0.8), width = 0.7, na.rm = TRUE) +
    geom_text(aes(label = ifelse(is.na(rappel), "0", sprintf("%.2f", rappel))),
              position = position_dodge(0.8), vjust = -0.4, size = 2.8, na.rm = TRUE) +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, 1) +
    labs(title = "Rappel par classe d'intensité (test temporel)",
         subtitle = "Part des cas réels de chaque état correctement prédits (1 = parfait)",
         x = NULL, y = "Rappel", fill = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  ggsave(file.path(outdir, "rappel_par_classe.png"), p, width = 8, height = 5, dpi = 120)
}


# =============================================================================
# POINT D'ENTRÉE : évaluation depuis le cache de covariables
# =============================================================================

#' Lancer l'évaluation complète à partir du panel en cache + config
run.model.evaluation <- function(config_path = "config.yaml") {
  cfg <- yaml::read_yaml(config_path)
  cache <- load.covariate.panel(cfg)
  if (is.null(cache)) stop("Aucun cache de covariables : lancez src/covariate_build.R")

  panel <- cache$data
  hex <- cache$hex
  cov_names <- names(cache$types)
  n_states <- length(levels(panel$target))

  # Paires de transition t -> t+1.
  pairs <- build.transition.pairs(panel, hex, cov_names)

  # Découpage temporel : par défaut, on réutilise celui du groupe 1 ; la section
  # `evaluation` du YAML peut le surcharger (phase de déclin = test_years).
  ev_cfg  <- cfg$evaluation %||% list()
  split   <- cfg$models$group1$temporal_split
  train_y <- ev_cfg$train_years %||% split$train_years
  test_y  <- ev_cfg$test_years  %||% split$test_years
  n_blocks <- ev_cfg$n_spatial_blocks %||% 4
  outdir  <- ev_cfg$outdir %||% "data/output/evaluation"

  # Libellés d'états issus des données réelles (Aucun/Léger/Modéré/Grave).
  labels <- levels(panel$target)

  spatiotemporal.evaluation(pairs, hex, cov_names, n_states,
                            train_years = train_y, test_years = test_y,
                            n_blocks = n_blocks, outdir = outdir, labels = labels,
                            cfg = cfg)
}

if (sys.nframe() == 0) {
  run.model.evaluation("config.yaml")
}
