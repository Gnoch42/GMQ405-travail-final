# =============================================================================
# src/covariate_evaluation.R
# -----------------------------------------------------------------------------
# OUTIL AUTONOME d'évaluation des covariables candidates.
#
# BUT : aider à CHOISIR les covariables AVANT de les intégrer au modèle. Cet
# outil est indépendant du pipeline final (mod_spa_temp.R ne l'appelle pas). On
# le lance ponctuellement, on lit les résultats, on décide quelles covariables
# garder, puis on les déclare dans config.yaml.
#
# Ce qu'il fait, pour une année de référence donnée :
#   1. Ramène la cible TBE et chaque covariable candidate sur la grille
#      hexagonale (réutilise src/features.R -> src/hex_grid.R).
#   2. Mesure l'ASSOCIATION de chaque covariable avec l'intensité TBE, avec le
#      test adapté au type de variable (continue / catégorielle / ordinale).
#   3. Détecte la COLINÉARITÉ entre covariables (redondance d'information).
#   4. Produit un TABLEAU classé + des GRAPHIQUES pour la décision.
#
# ---------------------------------------------------------------------------
# RAPPEL PÉDAGOGIQUE — les tests utilisés (pour non-spécialistes)
# ---------------------------------------------------------------------------
# * Corrélation de SPEARMAN (rho) : mesure si deux variables évoluent dans le
#   même sens, en se basant sur les RANGS (1er, 2e, 3e...) et non les valeurs
#   brutes. Vaut entre -1 et +1. |rho| proche de 1 = association forte. On
#   l'utilise pour une covariable continue (ou ordinale) vs l'intensité TBE
#   (elle-même ordinale), car elle ne suppose pas de relation linéaire stricte.
#
# * Test de KRUSKAL-WALLIS : version "sur les rangs" de l'ANOVA. Répond à :
#   « la covariable continue prend-elle des valeurs différentes selon la classe
#   d'intensité TBE ? ». Une p-value faible (< 0.05) suggère que oui (lien réel).
#
# * V de CRAMÉR : mesure l'association entre DEUX variables catégorielles
#   (ex. type de forêt vs classe d'intensité). Dérivé du test du khi-deux.
#   Vaut entre 0 (aucun lien) et 1 (lien parfait). Comparable au |rho| comme
#   force d'association, ce qui permet de classer ensemble covariables continues
#   et catégorielles.
#
# * COLINÉARITÉ / VIF : deux covariables colinéaires portent la même information
#   (ex. altitude et température très corrélées). En garder les deux nuit à
#   l'interprétation des modèles. Le VIF (Variance Inflation Factor) quantifie
#   cette redondance : VIF > 5 (voire 10) signale un problème. On calcule aussi
#   une matrice de corrélation entre covariables continues.
# =============================================================================

library(sf)
library(ggplot2)

# On réutilise la logique d'agrégation spatiale (tâche 2) et l'assemblage.
source("src/hex_grid.R")
source("src/features.R")


# -----------------------------------------------------------------------------
# MESURES D'ASSOCIATION ÉLÉMENTAIRES
# -----------------------------------------------------------------------------

#' V de Cramér entre deux variables catégorielles
#'
#' @return Liste : $strength (V dans [0,1]) et $p (p-value du khi-deux).
cramers.v <- function(x, y) {
  tab <- table(x, y)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(list(strength = NA_real_, p = NA_real_))
  chi <- suppressWarnings(chisq.test(tab))
  n <- sum(tab)
  v <- sqrt(as.numeric(chi$statistic) / (n * (min(dim(tab)) - 1)))
  list(strength = v, p = as.numeric(chi$p.value))
}


#' Associer UNE covariable à la cible ordinale
#'
#' Choisit automatiquement le test selon le type de covariable, et renvoie une
#' force d'association normalisée dans [0,1] pour permettre un classement commun.
#'
#' @param cov_vals   Valeurs de la covariable (une par hexagone).
#' @param target_rank Cible convertie en rang entier (1,2,3,...).
#' @param target_fac  Cible en facteur (pour Cramér's V).
#' @param value_type "continuous" | "ordinal" | "categorical".
#'
#' @return Liste : test, strength [0,1], p, direction.
associate.covariate <- function(cov_vals, target_rank, target_fac, value_type) {

  if (value_type %in% c("continuous", "ordinal")) {
    # Spearman : force = |rho| ; on garde aussi le signe (direction du lien).
    ok <- !is.na(cov_vals) & !is.na(target_rank)
    if (sum(ok) < 3) return(list(test = "spearman", strength = NA, p = NA, direction = NA))
    rho <- suppressWarnings(cor(as.numeric(cov_vals)[ok], target_rank[ok],
                                method = "spearman"))
    # Kruskal-Wallis pour la p-value (la covariable diffère-t-elle selon classe ?).
    kw <- suppressWarnings(tryCatch(
      kruskal.test(as.numeric(cov_vals)[ok], target_fac[ok])$p.value,
      error = function(e) NA_real_))
    list(test = "spearman", strength = abs(rho), p = kw,
         direction = ifelse(rho >= 0, "+", "-"))

  } else {
    # Catégorielle : V de Cramér vs classe d'intensité.
    cv <- cramers.v(cov_vals, target_fac)
    list(test = "cramers_v", strength = cv$strength, p = cv$p, direction = NA)
  }
}


# -----------------------------------------------------------------------------
# FONCTION PRINCIPALE
# -----------------------------------------------------------------------------

#' Évaluer un jeu de covariables candidates
#'
#' @param hex         Grille hexagonale (create.hex.grid()).
#' @param target_sf   sf de la cible TBE.
#' @param target_cfg  Config cible (year_field, value_field, ordinal_levels...).
#' @param covariates  Liste de covariables candidates (format YAML/simulé).
#' @param year        Année de référence pour l'analyse (coupe transversale).
#' @param min_coverage Seuil de couverture pour l'agrégation.
#' @param outdir      Dossier où écrire tableau (CSV) et graphiques (PNG).
#'
#' @return Liste : $ranking (data.frame classé), $collinearity (matrice),
#'   $vif (data.frame), $table (données hexagonales assemblées).
evaluate.covariates <- function(hex, target_sf, target_cfg, covariates,
                                year, min_coverage = 0.25,
                                outdir = "data/output/covariate_eval") {

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("== Évaluation des covariables (année de référence : ", year, ") ==")

  # --- 1. Agrégation cible + covariables sur la grille ---
  target_vals <- aggregate.target.year(hex, target_sf, target_cfg, year, min_coverage)
  # Absence d'infestation = niveau "Aucun" (information, pas donnée manquante).
  lv <- c("Aucun", target_cfg$ordinal_levels)
  target_fac <- factor(as.character(target_vals), levels = lv, ordered = TRUE)
  target_fac[is.na(target_fac)] <- "Aucun"
  target_rank <- as.integer(target_fac)   # 1=Aucun, 2=Léger, ...

  agg <- aggregate.covariates(hex, covariates, min_coverage)
  hex_tab <- st_drop_geometry(agg$hex)
  types <- agg$types
  cov_names <- names(covariates)

  # --- 2. Association de chaque covariable avec la cible ---
  rank_rows <- lapply(cov_names, function(nm) {
    a <- associate.covariate(hex_tab[[nm]], target_rank, target_fac, types[nm])
    data.frame(covariable = nm, type = types[nm], test = a$test,
               force = round(a$strength, 3), p_value = signif(a$p, 3),
               direction = ifelse(is.na(a$direction), "", a$direction),
               stringsAsFactors = FALSE)
  })
  ranking <- do.call(rbind, rank_rows)
  # Classement par force d'association décroissante (la plus prometteuse en tête).
  ranking <- ranking[order(-ranking$force), ]
  ranking$interpretation <- interpret.strength(ranking$force)

  # --- 3. Colinéarité entre covariables continues ---
  cont_names <- cov_names[types[cov_names] %in% c("continuous", "ordinal")]
  collinearity <- NULL
  vif_df <- NULL
  if (length(cont_names) >= 2) {
    cont_mat <- as.matrix(hex_tab[, cont_names, drop = FALSE])
    cont_mat <- apply(cont_mat, 2, as.numeric)
    collinearity <- cor(cont_mat, use = "pairwise.complete.obs", method = "spearman")
    vif_df <- compute.vif(as.data.frame(cont_mat))
  }

  # --- 4. Sorties : tableau CSV + graphiques PNG ---
  write.csv(ranking, file.path(outdir, "classement_covariables.csv"), row.names = FALSE)
  plot.ranking(ranking, outdir)
  if (!is.null(collinearity)) plot.collinearity(collinearity, outdir)

  message("Résultats écrits dans : ", outdir)
  print(ranking, row.names = FALSE)
  if (!is.null(vif_df)) { message("\n-- VIF (colinéarité, >5 = à surveiller) --"); print(vif_df, row.names = FALSE) }

  invisible(list(ranking = ranking, collinearity = collinearity,
                 vif = vif_df, table = hex_tab))
}


# -----------------------------------------------------------------------------
# UTILITAIRES DE SORTIE
# -----------------------------------------------------------------------------

#' Traduire une force d'association [0,1] en libellé lisible.
interpret.strength <- function(s) {
  cut(s, breaks = c(-Inf, 0.1, 0.3, 0.5, Inf),
      labels = c("négligeable", "faible", "modérée", "forte"))
}


#' VIF (Variance Inflation Factor) calculé "à la main"
#'
#' Pour chaque covariable continue, on régresse cette covariable sur toutes les
#' autres. Plus les autres l'expliquent (R² élevé), plus le VIF = 1/(1-R²) est
#' grand : c'est le signe qu'elle est redondante. Calcul sans package externe.
compute.vif <- function(df) {
  vars <- names(df)
  vif <- sapply(vars, function(v) {
    others <- setdiff(vars, v)
    fml <- as.formula(paste0("`", v, "` ~ ", paste0("`", others, "`", collapse = " + ")))
    r2 <- summary(lm(fml, data = df))$r.squared
    if (r2 >= 1) Inf else 1 / (1 - r2)
  })
  data.frame(covariable = vars, VIF = round(vif, 2), row.names = NULL)
}


#' Graphique en barres du classement des covariables.
plot.ranking <- function(ranking, outdir) {
  ranking$covariable <- factor(ranking$covariable, levels = rev(ranking$covariable))
  p <- ggplot(ranking, aes(x = covariable, y = force, fill = type)) +
    geom_col() +
    coord_flip() +
    labs(title = "Force d'association avec l'intensité TBE",
         subtitle = "Spearman |rho| (continu) ou V de Cramér (catégoriel), 0 = nul, 1 = fort",
         x = NULL, y = "Force d'association [0-1]") +
    theme_minimal(base_size = 12)
  ggsave(file.path(outdir, "classement_covariables.png"), p,
         width = 8, height = 5, dpi = 120)
}


#' Heatmap de colinéarité entre covariables continues (corrplot).
plot.collinearity <- function(cor_mat, outdir) {
  if (!requireNamespace("corrplot", quietly = TRUE)) return(invisible())
  png(file.path(outdir, "colinearite_covariables.png"), width = 700, height = 700)
  corrplot::corrplot(cor_mat, method = "color", type = "upper",
                     addCoef.col = "black", tl.col = "black",
                     title = "Colinéarité (Spearman) entre covariables continues",
                     mar = c(0, 0, 2, 0))
  dev.off()
}


# -----------------------------------------------------------------------------
# DÉMONSTRATION sur données SIMULÉES
# -----------------------------------------------------------------------------
# Permet de lancer l'outil immédiatement :  Rscript src/covariate_evaluation.R
# Génère un jeu factice, l'évalue, et écrit les sorties dans data/output/.
demo.covariate.evaluation <- function() {
  source("src/simulate_data.R")
  sim <- simulate.tbe.dataset(seed = 42, n_years = 8, extent_km = 200)

  hex <- create.hex.grid(sim$zone, cellsize = 15000, target_crs = 32198)

  target_cfg <- list(year_field = "year", value_field = "intensity",
                     value_type = "ordinal",
                     ordinal_levels = c("Léger", "Modéré", "Grave"))

  evaluate.covariates(
    hex = hex, target_sf = sim$tbe, target_cfg = target_cfg,
    covariates = sim$covariates, year = 2019, min_coverage = 0.25,
    outdir = "data/output/covariate_eval"
  )
}

# Exécution automatique uniquement si le fichier est lancé directement
# (Rscript src/covariate_evaluation.R), pas quand il est source() par un autre.
if (sys.nframe() == 0) {
  demo.covariate.evaluation()
}
