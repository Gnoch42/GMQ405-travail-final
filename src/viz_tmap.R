# =============================================================================
# src/viz_tmap.R
# -----------------------------------------------------------------------------
# Fonctions de CARTOGRAPHIE réutilisables (tmap v4), au STYLE UNIFORME.
#
# But : produire les cartes du projet (intensité TBE observée/prédite,
# covariables, gradients de risque) sans redéfinir le style à chaque fois. Une
# palette, une typographie et une mise en page de légende communes garantissent
# une esthétique professionnelle et cohérente d'une carte à l'autre.
#
# Deux fonctions principales, indépendantes du pipeline de modélisation (elles
# s'appliquent à N'IMPORTE quel objet sf du projet) :
#   - plot_donnees_brutes() : carte d'une variable telle quelle sur la grille
#                             hexagonale (continue ou catégorielle/ordinale).
#   - plot_flou_gaussien()  : lissage gaussien d'une variable -> SURFACE RASTER
#                             continue (gradients plus lisibles qu'un pavage brut).
#
# Dépendances : tmap (>= 4), sf, terra.
# =============================================================================

library(tmap)
library(sf)
library(terra)

`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# STYLE PARTAGÉ
# -----------------------------------------------------------------------------

#' Palette ORDONNÉE pour la sévérité TBE (4 états)
#'
#' Palette séquentielle de type "YlOrRd" (jaune -> rouge foncé, ColorBrewer),
#' standard professionnel lisible en impression et acceptable pour le daltonisme.
#' L'état "Absence" reçoit un gris neutre pour se distinguer nettement du gradient
#' d'intensité croissante.
PALETTE_TBE <- c("Absence" = "#f0f0f0", "Légère" = "#fed976",
                 "Modérée" = "#fd8d3c", "Sévère" = "#bd0026")

# Palette séquentielle par défaut pour les variables CONTINUES (viridis).
PALETTE_CONTINU <- "viridis"

#' Mise en page commune (typographie, légende, cadre) — thème du projet
#'
#' Regroupée ici pour que TOUTES les cartes partagent exactement le même style.
#' @return Un objet tmap (tm_layout) à ajouter aux cartes.
theme_tbe <- function() {
  tm_layout(
    text.fontfamily   = "sans",
    frame             = TRUE,
    frame.lwd         = 0.6,
    legend.outside    = TRUE,
    legend.outside.position = "right",
    legend.frame      = FALSE,
    legend.title.size = 0.95,
    legend.text.size  = 0.75,
    title.size        = 1.1,
    inner.margins     = c(0.02, 0.02, 0.02, 0.02)
  )
}

#' Sauvegarde optionnelle d'une carte (si `export` est fourni)
#' @param tm    Objet tmap.
#' @param export Chemin de sortie (png/pdf...) ou NULL (pas d'export).
#' @param dpi,width,height Options d'image (pouces).
export.carte <- function(tm, export, dpi = 300, width = 7, height = 6) {
  if (is.null(export)) return(invisible(NULL))
  dir.create(dirname(export), recursive = TRUE, showWarnings = FALSE)
  suppressMessages(tmap_save(tm, filename = export, dpi = dpi,
                             width = width, height = height, units = "in"))
  message("Carte enregistrée : ", export)
}


# -----------------------------------------------------------------------------
# 1. CARTE D'UNE VARIABLE BRUTE SUR LA GRILLE HEXAGONALE
# -----------------------------------------------------------------------------

#' Cartographier une variable telle quelle sur la grille hexagonale
#'
#' Détecte automatiquement le type de variable : catégorielle/ordinale (facteur
#' ou caractère) -> classes discrètes ; continue (numérique) -> gradient. La
#' variable "intensité TBE" (états Absence..Sévère) reçoit la palette dédiée.
#'
#' @param x        Objet sf (polygones hexagonaux).
#' @param var      Nom de la colonne à afficher.
#' @param titre    Titre de la carte (défaut : nom de la variable).
#' @param palette  Palette : vecteur de couleurs (catégoriel) ou nom de palette
#'                 (continu). NULL -> palette du projet selon le type.
#' @param style    Méthode de discrétisation du continu ("pretty", "quantile",
#'                 "jenks"...) transmise à tmap.
#' @param n        Nombre de classes pour le continu.
#' @param export   Chemin d'export (NULL = affichage seul).
#' @param dpi,width,height Options d'export.
#' @return L'objet tmap (invisible), affiché et/ou exporté.
plot_donnees_brutes <- function(x, var, titre = NULL, palette = NULL,
                                style = "pretty", n = 5, export = NULL,
                                dpi = 300, width = 7, height = 6) {
  stopifnot(inherits(x, "sf"), var %in% names(x))
  titre <- titre %||% var
  vals  <- x[[var]]
  is_cat <- is.factor(vals) || is.character(vals)

  if (is_cat) {
    # Catégoriel / ordinal : palette dédiée si les niveaux sont ceux de la TBE.
    lvls <- if (is.factor(vals)) levels(vals) else sort(unique(vals))
    pal  <- palette %||% if (all(lvls %in% names(PALETTE_TBE))) PALETTE_TBE[lvls] else "brewer.set2"
    scale <- if (is.ordered(vals)) tm_scale_ordinal(values = pal)
             else                  tm_scale_categorical(values = pal)
  } else {
    # Continu : gradient séquentiel discrétisé.
    pal   <- palette %||% PALETTE_CONTINU
    scale <- tm_scale_intervals(style = style, n = n, values = pal)
  }

  tm <- tm_shape(x) +
    tm_polygons(fill = var, fill.scale = scale,
                fill.legend = tm_legend(title = var),
                col = "grey40", lwd = 0.2) +
    tm_title(titre) +
    theme_tbe()

  export.carte(tm, export, dpi, width, height)
  print(tm)
  invisible(tm)
}


# -----------------------------------------------------------------------------
# 2. LISSAGE GAUSSIEN -> SURFACE RASTER CONTINUE
# -----------------------------------------------------------------------------

#' Lisser une variable et la cartographier en surface raster continue
#'
#' Convertit la grille hexagonale (pavage discret) en une surface continue plus
#' lisible : on rastérise la variable, puis on applique un flou gaussien
#' d'écart-type `sigma` (en unités de la carte, p. ex. mètres). Utile pour
#' visualiser des GRADIENTS (pression d'infestation, risque) sans l'effet
#' "mosaïque" des hexagones.
#'
#' @param x          Objet sf (hexagones) contenant une variable NUMÉRIQUE.
#' @param var        Nom de la variable continue à lisser.
#' @param sigma      Intensité du flou = écart-type du noyau gaussien, en unités
#'                   du CRS (mètres si CRS projeté). Plus grand = plus lisse.
#' @param resolution Taille de cellule du raster de sortie (mêmes unités).
#' @param titre,palette,export,dpi,width,height  Comme plot_donnees_brutes().
#' @return L'objet tmap (invisible). Le raster lissé est aussi renvoyé en attribut
#'   "raster" pour réutilisation éventuelle.
plot_flou_gaussien <- function(x, var, sigma, resolution = 250, titre = NULL,
                               palette = NULL, export = NULL,
                               dpi = 300, width = 7, height = 6) {
  stopifnot(inherits(x, "sf"), var %in% names(x))
  if (!is.numeric(x[[var]]))
    stop("plot_flou_gaussien() attend une variable NUMÉRIQUE (continue).")
  titre <- titre %||% paste0(var, " (lissé, sigma = ", sigma, ")")

  # 1. Gabarit raster couvrant l'emprise des hexagones, à la résolution demandée.
  template <- terra::rast(terra::ext(terra::vect(x)), resolution = resolution,
                          crs = terra::crs(terra::vect(x)))
  # 2. Rastérisation de la variable (valeur de l'hexagone -> cellules couvertes).
  r <- terra::rasterize(terra::vect(x), template, field = var)

  # 3. Noyau gaussien d'écart-type sigma (converti en cellules par focalMat).
  w <- terra::focalMat(r, d = sigma, type = "Gauss")
  smooth <- terra::focal(r, w = w, fun = "sum", na.rm = TRUE, na.policy = "omit")
  names(smooth) <- var

  # 4. Carte raster continue au style commun.
  pal <- palette %||% PALETTE_CONTINU
  tm <- tm_shape(smooth) +
    tm_raster(col.scale = tm_scale_continuous(values = pal),
              col.legend = tm_legend(title = var)) +
    tm_title(titre) +
    theme_tbe()

  export.carte(tm, export, dpi, width, height)
  print(tm)
  attr(tm, "raster") <- smooth
  invisible(tm)
}


# -----------------------------------------------------------------------------
# DÉMONSTRATION (Rscript src/viz_tmap.R) : cartes depuis le cache de covariables
# -----------------------------------------------------------------------------
demo.viz <- function(config_path = "config.yaml") {
  source("src/features.R")
  cfg <- yaml::read_yaml(config_path)
  cache <- load.covariate.panel(cfg)
  if (is.null(cache)) stop("Cache absent : lancez src/covariate_build.R d'abord.")

  hex <- cache$hex
  # On joint une covariable statique (prop_hote) à la géométrie pour l'exemple.
  stat <- unique(cache$data[, c("hex_id", "prop_hote")])
  hex <- merge(hex, stat, by = "hex_id", all.x = TRUE)

  tmap_mode("plot")
  plot_donnees_brutes(hex, "prop_hote", titre = "Proportion d'essences hôtes",
                      export = "data/output/viz/prop_hote_brut.png")
  plot_flou_gaussien(hex, "prop_hote", sigma = 3000, resolution = 250,
                     titre = "Proportion d'essences hôtes (lissée)",
                     export = "data/output/viz/prop_hote_flou.png")
}

if (sys.nframe() == 0) {
  demo.viz("config.yaml")
}
