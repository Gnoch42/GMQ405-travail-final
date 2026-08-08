# =============================================================================
# src/models_group2.R
# -----------------------------------------------------------------------------
# GROUPE 2 — Modèles économétriques/statistiques à UN INSTANT t.
#
# But de ce groupe : COMPRENDRE (inférence) les relations entre les covariables
# et l'intensité TBE à un moment donné. Il ne s'agit PAS de prédire le futur
# (c'est le rôle du groupe 1), mais d'estimer « quelle covariable est associée à
# quelle intensité, et l'espace joue-t-il un rôle ? ».
#
# Analyse en COUPE TRANSVERSALE (cross-section) : on fige une année, une ligne
# par hexagone. La cible ordinale est traitée ici comme numérique (rang 0,1,2,3)
# pour rester dans le cadre des modèles linéaires spatiaux classiques — choix
# assumé et documenté ; une alternative serait un modèle ordinal (MASS::polr),
# laissé en extension.
#
# Les 7 modèles, du plus simple au plus riche :
#   MCO           : régression linéaire ordinaire, IGNORE l'espace (référence).
#   SLX           : MCO + moyennes des covariables chez les VOISINS (spillover X).
#   SAR (lag)     : la valeur d'un hexagone dépend de la valeur de ses VOISINS.
#   SEM (erreur)  : l'espace agit via des ERREURS corrélées entre voisins.
#   Durbin (SDM)  : SAR + covariables des voisins (le plus général en "lag").
#   Durbin erreur : SEM + covariables des voisins (le plus général en "erreur").
#   GAM           : relations NON linéaires lissées (splines) entre X et cible.
#
# ---------------------------------------------------------------------------
# GLOSSAIRE (pour non-spécialistes)
# ---------------------------------------------------------------------------
# * MATRICE DE VOISINAGE / CONTIGUÏTÉ : tableau qui dit quels hexagones sont
#   "voisins". Ici, deux hexagones sont voisins s'ils se touchent. On la
#   normalise pour que chaque hexagone reçoive une MOYENNE de ses voisins.
# * AUTOCORRÉLATION SPATIALE : tendance de lieux proches à se ressembler. Le
#   test de MORAN mesure cela sur les résidus : s'il reste de l'autocorrélation
#   dans les résidus du MCO, c'est que l'espace n'est pas correctement modélisé
#   (justifie les modèles spatiaux SAR/SEM/Durbin).
# * AIC / BIC : critères comparant les modèles ; plus BAS = meilleur compromis
#   ajustement/complexité. Comparables ENTRE modèles du même groupe.
#
# Dépendances : spdep, spatialreg, mgcv (tous installés).
# =============================================================================

library(sf)


#' Construire la matrice de voisinage (contiguïté) des hexagones
#'
#' @param hex Grille hexagonale (avec géométrie).
#' @return Un objet `listw` (spdep) : voisinage normalisé en ligne (chaque
#'   hexagone pondère ses voisins pour en faire une moyenne). style = "W".
build.spatial.weights <- function(hex) {
  # poly2nb : détecte les voisins par contiguïté (hexagones qui se touchent).
  nb <- spdep::poly2nb(hex, queen = TRUE)
  # nb2listw : transforme en poids ; zero.policy autorise les hexagones isolés.
  spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
}


#' p-value du test de Moran sur les RÉSIDUS d'un modèle (tout type)
#'
#' Diagnostic uniforme appliqué aux sept modèles : reste-t-il de l'autocorrélation
#' spatiale dans les résidus ? Pour le MCO, une valeur très faible justifie le
#' passage aux modèles spatiaux ; pour les modèles spatiaux (SAR/SEM/Durbin), une
#' valeur REDEVENUE ÉLEVÉE (> 0,05) confirme qu'ils ont bien absorbé la structure
#' spatiale. On applique le test de Moran directement aux résidus (`moran.test`),
#' ce qui fonctionne pour lm, les objets `spatialreg` et `mgcv`.
#'
#' @param model Modèle ajusté (lm, Sarlm, SLX, gam...).
#' @param listw Matrice de voisinage (listw) alignée sur les résidus.
#' @return La p-value (NA si le calcul échoue).
residual.moran.p <- function(model, listw) {
  tryCatch(
    spdep::moran.test(as.numeric(residuals(model)), listw,
                      zero.policy = TRUE, na.action = stats::na.omit)$p.value,
    error = function(e) NA_real_)
}


#' Ajuster les 7 modèles du groupe 2 pour une année donnée
#'
#' @param data    data.frame d'UNE année : colonnes hex_id, target (ordinale),
#'                + covariables. Produit par build.hex.panel() filtré sur l'année.
#' @param hex     Grille hexagonale (géométrie) alignée sur `data` (même hex_id).
#' @param cov_names Noms des covariables à utiliser comme régresseurs.
#' @param models  Sous-ensemble de modèles à ajuster (défaut : les 7).
#'
#' @return data.frame standardisé : une ligne par modèle, colonnes de métriques
#'   comparables (voir standardize.group2.result()).
fit.group2 <- function(data, hex, cov_names,
                       models = c("ols","slx","sar","sem","durbin","durbin_error","gam")) {

  # Garde-fou : les modèles de régression exigent au moins une covariable.
  # Sans covariable déclarée dans le YAML, on ne peut rien estimer -> on renvoie
  # une ligne d'information plutôt que de planter le pipeline.
  if (length(cov_names) == 0) {
    message("  [Groupe 2] aucune covariable déclarée -> groupe ignoré ",
            "(ajoutez des covariables dans config `covariates:`).")
    return(standardize.group2.result("(aucune covariable)", NULL, statut = "aucune_covariable"))
  }

  # --- 1. Préparation : cible ordinale -> numérique (rang) ---
  # Nécessaire pour les modèles linéaires spatiaux. On garde l'ordre Aucun < ...
  data$y <- as.integer(data$target) - 1L   # 0 = Aucun, 1 = Léger, ...

  # On ne garde que les hexagones complets (cible + toutes covariables non NA),
  # sinon les modèles spatiaux (qui exigent une matrice de voisinage cohérente)
  # échouent. On aligne la grille sur ces mêmes hexagones.
  keep <- stats::complete.cases(data[, c("y", cov_names)])
  data <- data[keep, ]
  hex  <- hex[match(data$hex_id, hex$hex_id), ]

  listw <- build.spatial.weights(hex)

  # Formule de base : y ~ cov1 + cov2 + ...
  fml <- as.formula(paste("y ~", paste(cov_names, collapse = " + ")))

  results <- list()

  # --- MCO (référence non spatiale) ---
  if ("ols" %in% models) {
    m <- lm(fml, data = data)
    results$ols <- standardize.group2.result("MCO", m,
      aic = AIC(m), bic = BIC(m), r2 = summary(m)$r.squared,
      moran_p = residual.moran.p(m, listw), loglik = as.numeric(logLik(m)))
  }

  # --- SLX : ajoute les covariables moyennes des voisins (lagged X) ---
  if ("slx" %in% models) {
    m <- tryCatch(spatialreg::lmSLX(fml, data = data, listw = listw, zero.policy = TRUE),
                  error = function(e) NULL)
    results$slx <- safe.result("SLX", m, listw)
  }

  # --- SAR (lag) : y dépend de la moyenne de y chez les voisins ---
  if ("sar" %in% models) {
    m <- tryCatch(spatialreg::lagsarlm(fml, data = data, listw = listw, zero.policy = TRUE),
                  error = function(e) NULL)
    results$sar <- safe.result("SAR (lag)", m, listw)
  }

  # --- SEM (erreur) : l'espace agit via des erreurs corrélées ---
  if ("sem" %in% models) {
    m <- tryCatch(spatialreg::errorsarlm(fml, data = data, listw = listw, zero.policy = TRUE),
                  error = function(e) NULL)
    results$sem <- safe.result("SEM (erreur)", m, listw)
  }

  # --- Durbin spatial (SDM) : SAR + X des voisins ---
  if ("durbin" %in% models) {
    m <- tryCatch(spatialreg::lagsarlm(fml, data = data, listw = listw,
                                       Durbin = TRUE, zero.policy = TRUE),
                  error = function(e) NULL)
    results$durbin <- safe.result("Durbin spatial", m, listw)
  }

  # --- Durbin erreur (SDEM) : SEM + X des voisins ---
  if ("durbin_error" %in% models) {
    m <- tryCatch(spatialreg::errorsarlm(fml, data = data, listw = listw,
                                         Durbin = TRUE, zero.policy = TRUE),
                  error = function(e) NULL)
    results$durbin_error <- safe.result("Durbin erreur", m, listw)
  }

  # --- GAM : relations non linéaires (splines) sur chaque covariable continue ---
  if ("gam" %in% models) {
    # s(x) = spline lissée ; on l'applique aux covariables numériques, les
    # catégorielles restent en effet linéaire classique.
    is_num <- sapply(cov_names, function(v) is.numeric(data[[v]]))
    terms <- ifelse(is_num, paste0("s(", cov_names, ")"), cov_names)
    gfml <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
    m <- tryCatch(mgcv::gam(gfml, data = data), error = function(e) NULL)
    if (!is.null(m)) {
      results$gam <- standardize.group2.result("GAM", m,
        aic = AIC(m), bic = BIC(m), r2 = summary(m)$r.sq,
        moran_p = residual.moran.p(m, listw), loglik = as.numeric(logLik(m)))
    } else {
      results$gam <- standardize.group2.result("GAM", NULL)
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


# -----------------------------------------------------------------------------
# STANDARDISATION DE LA SORTIE (structure commune au GROUPE 2)
# -----------------------------------------------------------------------------
# Toutes les métriques ci-dessous sont COMPARABLES entre les 7 modèles :
#   - AIC / BIC : plus bas = meilleur (compromis ajustement/complexité).
#   - pseudo-R2 : part de variance expliquée (indicatif).
#   - moran_p   : p-value de Moran sur les résidus. Une valeur ÉLEVÉE (> 0.05)
#                 est BONNE ici : elle indique qu'il ne reste plus
#                 d'autocorrélation spatiale non modélisée.
#   - rho/lambda : coefficients spatiaux (force de la dépendance spatiale).

standardize.group2.result <- function(name, model,
                                       aic = NA_real_, bic = NA_real_,
                                       r2 = NA_real_, moran_p = NA_real_,
                                       rho = NA_real_, lambda = NA_real_,
                                       loglik = NA_real_, statut = "ok") {
  # Un modèle NULL signale un échec d'ajustement, sauf si l'appelant a déjà
  # fourni un statut plus précis (ex. "aucune_covariable").
  if (is.null(model) && statut == "ok") statut <- "échec_ajustement"
  data.frame(
    modele    = name,
    statut    = statut,
    AIC       = round(aic, 1),
    BIC       = round(bic, 1),
    pseudo_R2 = round(r2, 3),
    moran_p_residus = signif(moran_p, 3),
    rho       = round(rho, 3),
    lambda    = round(lambda, 3),
    logLik    = round(loglik, 1),
    stringsAsFactors = FALSE
  )
}


#' Extraire les métriques standard d'un modèle spatialreg (SAR/SEM/Durbin/SLX)
safe.result <- function(name, model, listw) {
  if (is.null(model)) return(standardize.group2.result(name, NULL))

  s <- tryCatch(summary(model), error = function(e) NULL)
  aic <- tryCatch(AIC(model), error = function(e) NA_real_)
  bic <- tryCatch(BIC(model), error = function(e) NA_real_)
  ll  <- tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_)
  # rho (dépendance sur y) et lambda (sur l'erreur) selon le type de modèle.
  rho    <- if (!is.null(model$rho)) model$rho else NA_real_
  lambda <- if (!is.null(model$lambda)) model$lambda else NA_real_
  # pseudo-R2 : corrélation² entre valeurs ajustées et observées.
  # suppressMessages : spatialreg émet un avis "response is known" sans incidence.
  r2 <- suppressMessages(tryCatch({
    fitted_vals <- as.numeric(fitted(model))
    cor(fitted_vals, model$y)^2
  }, error = function(e) NA_real_))

  standardize.group2.result(name, model, aic = aic, bic = bic, r2 = r2,
                            moran_p = residual.moran.p(model, listw),
                            rho = rho, lambda = lambda, loglik = ll)
}
