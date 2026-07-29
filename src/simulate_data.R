# =============================================================================
# src/simulate_data.R
# -----------------------------------------------------------------------------
# Générateur de DONNÉES FACTICES pour la démonstration du pipeline.
#
# But : permettre de lancer et tester TOUT le pipeline (grille, agrégation,
# évaluation des covariables, modélisation) SANS attendre les vraies données.
# Dès que les chemins réels sont renseignés dans config.yaml, ces fonctions ne
# sont plus appelées (voir le drapeau `simulate$use_simulated`).
#
# Ce que l'on simule :
#   - Une zone d'étude (contour rectangulaire simple, en Lambert québécois).
#   - La cible TBE : des polygones d'intensité ordinale (Léger/Modéré/Grave),
#     avec un FOYER qui se déplace et grossit d'année en année (pour imiter
#     une dynamique de propagation spatio-temporelle).
#   - Des covariables factices :
#       * un RASTER continu (ex. "température", corrélé à la latitude) ;
#       * un RASTER continu bruité (covariable peu informative, pour tester la
#         détection de non-pertinence) ;
#       * des POLYGONES catégoriels (ex. "type de forêt").
#
# Toutes les géométries sont produites directement dans un CRS projeté en
# mètres (EPSG:32198, Lambert québécois) pour rester cohérent avec la grille.
# =============================================================================

library(sf)
library(terra)

# CRS de travail commun aux données simulées (projeté, en mètres).
.SIM_CRS <- 32198


#' Générer un jeu de données TBE + covariables simulé
#'
#' @param seed    Graine aléatoire (reproductibilité).
#' @param n_years Nombre d'années à simuler (le foyer progresse chaque année).
#' @param extent_km Côté de la zone d'étude carrée, en kilomètres.
#'
#' @return Une liste avec :
#'   - $zone       : sf, contour de la zone d'étude.
#'   - $tbe        : sf, polygones TBE (colonnes `year`, `intensity` ordinale).
#'   - $covariates : liste nommée de sources de covariables, chacune une liste
#'                   {source, type, value_type, value_field} directement
#'                   consommable par aggregate.to.hex().
simulate.tbe.dataset <- function(seed = 42, n_years = 12, extent_km = 200) {
  set.seed(seed)

  # --- 1. Zone d'étude : un carré de `extent_km` de côté ---
  # Coordonnées arbitraires plausibles en Lambert québécois (mètres).
  x0 <- 100000; y0 <- 200000
  side <- extent_km * 1000
  zone <- st_sf(
    geometry = st_sfc(st_polygon(list(rbind(
      c(x0, y0), c(x0 + side, y0),
      c(x0 + side, y0 + side), c(x0, y0 + side), c(x0, y0)
    ))), crs = .SIM_CRS)
  )

  # --- 2. Cible TBE : un foyer qui se déplace et grossit chaque année ---
  # Modèle simple : un centre de foyer part du coin sud-ouest et migre vers le
  # nord-est ; son rayon augmente. Trois anneaux concentriques donnent les trois
  # niveaux d'intensité (cœur = Grave, périphérie = Léger).
  levels_ord <- c("Léger", "Modéré", "Grave")
  years <- seq_len(n_years) + 2013            # 2014, 2015, ...
  tbe_list <- list()

  for (k in seq_along(years)) {
    frac <- (k - 1) / max(1, n_years - 1)     # progression 0 -> 1
    cx <- x0 + side * (0.25 + 0.5 * frac)     # déplacement du foyer
    cy <- y0 + side * (0.25 + 0.5 * frac)
    base_r <- side * (0.08 + 0.12 * frac)     # rayon croissant

    centre <- st_sfc(st_point(c(cx, cy)), crs = .SIM_CRS)
    # Trois anneaux : Grave (petit), Modéré (moyen), Léger (grand).
    ring_grave  <- st_buffer(centre, base_r * 0.5)
    ring_modere <- st_difference(st_buffer(centre, base_r * 1.0), ring_grave)
    ring_leger  <- st_difference(st_buffer(centre, base_r * 1.6),
                                 st_buffer(centre, base_r * 1.0))

    yr_sf <- st_sf(
      year      = years[k],
      intensity = factor(c("Grave", "Modéré", "Léger"), levels = levels_ord),
      geometry  = c(ring_grave, ring_modere, ring_leger)
    )
    # On garde uniquement la partie dans la zone d'étude.
    yr_sf <- suppressWarnings(st_intersection(yr_sf, st_geometry(zone)))
    tbe_list[[k]] <- yr_sf
  }
  tbe <- do.call(rbind, tbe_list)
  st_geometry(tbe) <- "geometry"

  # --- 3. Covariables factices ---
  covariates <- list()

  # 3a. Raster "température" : gradient nord-sud (corrélé à la position, donc
  #     partiellement corrélé au déplacement du foyer -> covariable pertinente).
  r_template <- rast(ext(x0, x0 + side, y0, y0 + side),
                     resolution = side / 100, crs = paste0("EPSG:", .SIM_CRS))
  yy <- yFromCell(r_template, seq_len(ncell(r_template)))
  temp_vals <- 20 - 15 * (yy - y0) / side + rnorm(ncell(r_template), 0, 1)
  temp <- setValues(r_template, temp_vals)
  names(temp) <- "temperature"
  covariates[["temperature"]] <- list(
    source = temp, type = "raster", value_type = "continuous", field = NULL)

  # 3b. Raster "bruit" : pur bruit blanc, AUCUN lien avec la cible.
  #     Sert à vérifier que l'outil d'évaluation le classe bien comme non pertinent.
  noise <- setValues(r_template, rnorm(ncell(r_template)))
  names(noise) <- "bruit"
  covariates[["bruit"]] <- list(
    source = noise, type = "raster", value_type = "continuous", field = NULL)

  # 3c. Polygones "type de forêt" : découpage aléatoire de la zone en cellules
  #     carrées, chacune recevant une catégorie. Corrélation modérée avec la
  #     latitude (donc lien indirect avec la cible).
  grid_poly <- st_make_grid(zone, cellsize = side / 12, square = TRUE)
  grid_poly <- st_sf(geometry = grid_poly)
  cy_poly <- st_coordinates(st_centroid(grid_poly))[, 2]
  prob_conifere <- (cy_poly - y0) / side          # + on va au nord, + de conifères
  grid_poly$foret <- ifelse(runif(nrow(grid_poly)) < prob_conifere,
                            "Conifere", "Feuillu")
  grid_poly$foret <- factor(grid_poly$foret)
  covariates[["type_foret"]] <- list(
    source = grid_poly, type = "vector", value_type = "categorical",
    field = "foret")

  list(zone = zone, tbe = tbe, covariates = covariates)
}
