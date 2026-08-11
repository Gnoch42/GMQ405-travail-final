# =============================================================================
# src/visualisation.R
# -----------------------------------------------------------------------------
# FONCTIONS GRAPHIQUES DU PROJET — point unique pour TOUS les visuels.
#
# Ce fichier regroupe l'ensemble des fonctions de cartographie (tmap) et de
# graphiques (ggplot2) du projet, afin de garantir une esthétique UNIFORME :
# une palette d'intensité TBE commune, une palette de typologie LISA commune, et
# deux thèmes partagés (un pour les cartes tmap, un pour les graphiques ggplot).
#
# Organisation :
#   0. Palettes et thèmes partagés
#   1. Cartes génériques d'une variable (brute / lissée) sur la grille hexagonale
#   2. Cartes de l'ampleur de l'épidémie (intensité observée, animation)
#   3. Cartes de dynamique spatiale (typologie LISA)
#   4. Graphiques d'aide à la décision (choix de matrice, optimisation geocmeans)
#   5. Graphiques d'évaluation des covariables (classement, colinéarité)
#   6. Graphiques d'évaluation des modèles (confusion, kappa, rappel par classe)
#
# Convention : chaque fonction qui enregistre un fichier accepte un dossier de
# sortie (`outdir`) ; les valeurs par défaut pointent vers le dossier `outputs/`
# subdivisé par étape.
#
# Dépendances : tmap (>= 4), sf, terra, ggplot2, ggpubr, tidyr (+ corrplot, gifski).
# =============================================================================

library(tmap)
library(sf)
library(terra)
library(ggplot2)
library(ggpubr)
library(tidyr)

`%||%` <- function(a, b) if (is.null(a)) b else a


# 0. PALETTES ET THÈMES PARTAGÉS -------------------------------------------

#' Palette ORDONNÉE de l'intensité TBE.
#'
#' Gradient séquentiel jaune -> rouge foncé (type ColorBrewer "YlOrRd"), lisible
#' en impression et acceptable pour le daltonisme ; l'absence reçoit un gris
#' neutre. Le vecteur couvre les DEUX conventions de libellés utilisées dans le
#' projet ("Aucun/Léger/Modéré/Grave" et "Absence/Légère/Modérée/Sévère"), de
#' sorte que n'importe quelle carte d'intensité utilise les mêmes couleurs.
PALETTE_TBE <- c(
  "Aucun"   = "#f0f0f0", "Absence" = "#f0f0f0",
  "Léger"   = "#ffffb2", "Légère"  = "#ffffb2",
  "Modéré"  = "#fd8d3c", "Modérée" = "#fd8d3c",
  "Grave"   = "#bd0026", "Sévère"  = "#bd0026"
)

# Palette séquentielle par défaut pour les variables CONTINUES (viridis).
PALETTE_CONTINU <- "viridis"

# Libellés des quatre états ordonnés (index 0..3) — partagés par l'évaluation.
ETATS_TBE <- c("Aucun", "Léger", "Modéré", "Grave")

#' Palette des typologies d'autocorrélation locale (nuage de Moran / LISA).
PALETTE_TYPO_LISA <- c("HH" = "red", "HL" = "lightpink", "LL" = "blue",
                       "LH" = "skyblue2", "Non sign." = "lightgray")

#' Thème commun des CARTES (tmap) : typographie, légende, cadre.
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

#' Thème commun des GRAPHIQUES (ggplot2) : minimal, légende en haut.
theme_tbe_gg <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey30"))
}

#' Enregistrer une carte tmap si un chemin `export` est fourni (sinon rien).
export.carte <- function(tm, export, dpi = 300, width = 7, height = 6) {
  if (is.null(export)) return(invisible(NULL))
  dir.create(dirname(export), recursive = TRUE, showWarnings = FALSE)
  suppressMessages(tmap_save(tm, filename = export, dpi = dpi,
                             width = width, height = height, units = "in"))
  message("Carte enregistrée : ", export)
}

#' Enregistrer un graphique ggplot si un chemin `export` est fourni (sinon rien).
export.graphique <- function(p, export, dpi = 120, width = 8, height = 5) {
  if (is.null(export)) return(invisible(NULL))
  dir.create(dirname(export), recursive = TRUE, showWarnings = FALSE)
  ggsave(export, p, width = width, height = height, dpi = dpi)
  message("Graphique enregistré : ", export)
}


# 1. CARTES GÉNÉRIQUES D'UNE VARIABLE SUR LA GRILLE HEXAGONALE -------------

#' Cartographier une variable telle quelle sur la grille hexagonale
#'
#' Détecte le type de variable : facteur ORDONNÉ (intensité) -> palette
#' séquentielle ; facteur/caractère non ordonné -> palette catégorielle ;
#' numérique -> gradient continu discrétisé.
#'
#' @param x        Objet sf (polygones hexagonaux).
#' @param var      Nom de la colonne à afficher.
#' @param titre    Titre (défaut : nom de la variable).
#' @param palette  Palette explicite (prioritaire) : vecteur de couleurs ou nom.
#' @param style,n  Discrétisation du continu (méthode et nombre de classes).
#' @param export   Chemin d'export (NULL = affichage seul).
#' @param dpi,width,height Options d'export.
#' @return L'objet tmap (invisible).
plot_donnees_brutes <- function(x, var, titre = NULL, palette = NULL,
                                style = "pretty", n = 5, export = NULL,
                                dpi = 300, width = 7, height = 6) {
  stopifnot(inherits(x, "sf"), var %in% names(x))
  titre <- titre %||% var
  vals  <- x[[var]]
  is_cat <- is.factor(vals) || is.character(vals)

  if (is_cat) {
    lvls <- if (is.factor(vals)) levels(vals) else sort(unique(vals))
    if (is.ordered(vals)) {
      # Ordinal (intensité) : palette séquentielle, indépendante des libellés.
      pal <- palette %||% "brewer.yl_or_rd"
      scale <- tm_scale_ordinal(values = pal)
    } else {
      # Catégoriel non ordonné : palette TBE si applicable, sinon Set2.
      pal <- palette %||% if (all(lvls %in% names(PALETTE_TBE))) unname(PALETTE_TBE[lvls]) else "brewer.set2"
      scale <- tm_scale_categorical(values = pal)
    }
  } else {
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


#' Lisser une variable et la cartographier en surface raster continue
#'
#' Rastérise la variable puis applique un flou gaussien d'écart-type `sigma`
#' (unités du CRS). Utile pour visualiser des GRADIENTS sans l'effet "mosaïque"
#' des hexagones.
#'
#' @param x          Objet sf (hexagones) avec une variable NUMÉRIQUE.
#' @param var        Nom de la variable continue à lisser.
#' @param sigma      Écart-type du noyau gaussien (unités du CRS ; + grand = + lisse).
#' @param resolution Taille de cellule du raster de sortie.
#' @param titre,palette,export,dpi,width,height  Comme plot_donnees_brutes().
#' @return L'objet tmap (invisible) ; le raster lissé est en attribut "raster".
plot_flou_gaussien <- function(x, var, sigma, resolution = 250, titre = NULL,
                               palette = NULL, export = NULL,
                               dpi = 300, width = 7, height = 6) {
  stopifnot(inherits(x, "sf"), var %in% names(x))
  if (!is.numeric(x[[var]]))
    stop("plot_flou_gaussien() attend une variable NUMÉRIQUE (continue).")
  titre <- titre %||% paste0(var, " (lissé, sigma = ", sigma, ")")

  template <- terra::rast(terra::ext(terra::vect(x)), resolution = resolution,
                          crs = terra::crs(terra::vect(x)))
  r <- terra::rasterize(terra::vect(x), template, field = var)
  w <- terra::focalMat(r, d = sigma, type = "Gauss")
  smooth <- terra::focal(r, w = w, fun = "sum", na.rm = TRUE, na.policy = "omit")
  names(smooth) <- var

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


# 2. CARTES DE L'AMPLEUR DE L'ÉPIDÉMIE -------------------------------------

#' Carte de l'intensité de défoliation observée (polygones TBE sur la zone)
#'
#' @param tbe        sf des polygones TBE (colonne "Niveau" ordinale).
#' @param study.zone sf de la zone d'étude (fond de carte).
#' @return L'objet tmap.
map.tbe <- function(tbe, study.zone) {
  tmap_mode("plot")
  tm_shape(study.zone) +
    tm_polygons(fill = "grey95", col = "black", lwd = 1) +
    tm_shape(tbe) +
    tm_polygons(
      col = NA, fill = "Niveau",
      fill.scale  = tm_scale_categorical(values = PALETTE_TBE[levels(tbe$Niveau)]),
      fill.legend = tm_legend(title = "Intensité de la défoliation",
                              frame = FALSE, position = c(0.7, 0.2))
    ) +
    tm_layout(frame = FALSE) +
    tm_scalebar(position = tm_pos_in("right", "bottom"),
                breaks = c(0, 20), text.size = 0.7) +
    tm_credits(
      "Source : MRNF, Données Québec\nAuteur : Jean Yves, Antoine Fortier, Rémy Billette",
      position = tm_pos_out("center", "bottom")
    )
}

#' Animation de l'évolution annuelle de l'épidémie (GIF)
#'
#' @param outdir Dossier d'écriture du GIF (défaut : outputs/ampleur).
map.tbe.animation <- function(tbe, study.zone, outdir = "outputs/ampleur") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  tbe.clip <- st_intersection(tbe, study.zone)
  carte <- tm_shape(study.zone) +
    tm_polygons(fill = "grey95", col = "black") +
    tm_shape(tbe.clip) +
    tm_polygons(
      col = NA, fill = "Niveau",
      fill.scale  = tm_scale_categorical(values = PALETTE_TBE[levels(tbe.clip$Niveau)]),
      fill.legend = tm_legend(title = "Intensité", frame = FALSE)
    ) +
    tm_animate(frames = "ANNEE", fps = 1L, play = "loop") +
    tm_title("Évolution de l'épidémie de TBE") +
    tm_layout(frame = FALSE)
  tmap_animation(carte, filename = file.path(outdir, "animation_tbe.gif"))
  carte
}
# Alias rétro-compatible (ancien nom).
tbe.animation <- function(tbe, study.zone, ...) map.tbe.animation(tbe, study.zone, ...)


# 3. CARTES DE DYNAMIQUE SPATIALE (TYPOLOGIE LISA) ---------------------------

#' Carte de la typologie LISA (avec ou sans légende)
#'
#' @param data   sf contenant le champ de typologie (HH/HL/LL/LH/Non sign.).
#' @param field  Nom de la colonne de typologie.
#' @param year   Année (utilisée pour le titre et le nom de fichier).
#' @param legend TRUE = avec légende ; FALSE = sans (pour les planches multi-années).
#' @param outdir Dossier d'écriture (défaut : outputs/dynamique).
lisa.map <- function(data, field, year, legend = TRUE, outdir = "outputs/dynamique") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  leg <- if (legend) tm_legend(title = "Typologie",
                               position = tm_pos_in("right", "bottom"), frame = FALSE)
         else tm_legend_hide()
  map <- tm_shape(data) +
    tm_polygons(col = "black", lwd = 0.5, fill = field,
                fill.scale = tm_scale_categorical(values = PALETTE_TYPO_LISA),
                fill.legend = leg) +
    tm_title(as.character(year)) +
    tm_layout(frame = FALSE)
  if (legend) {
    map <- map + tm_scalebar(breaks = c(0, 20), text.size = 1.2,
                             position = tm_pos_in("RIGHT", "BOTTOM"))
  }
  suffix <- if (legend) "legend" else ""
  tmap_save(map, file.path(outdir, paste0("lisa_", year, suffix, ".png")),
            width = 5, height = 4, dpi = 300)
  map
}
# Alias rétro-compatible (ancien nom sans légende).
lisa.map.no.legend <- function(data, field, year, outdir = "outputs/dynamique") {
  lisa.map(data, field, year, legend = FALSE, outdir = outdir)
}


# 4. GRAPHIQUES D'AIDE À LA DÉCISION ---------------------------------------

#' Comparaison des matrices de pondération spatiale (I de Moran)
moran.I.comparison <- function(df, export = NULL) {
  p <- ggplot(df, aes(x = reorder(mat, moran.I), y = moran.I)) +
    geom_segment(aes(xend = reorder(mat, moran.I), y = 0, yend = moran.I)) +
    geom_point(size = 4, fill = "red", shape = 21) +
    labs(x = "Matrice de pondération spatiale", y = "I de Moran") +
    coord_flip() + theme_tbe_gg()
  export.graphique(p, export); p
}

#' Optimisation du paramètre m (geocmeans) : inertie et silhouette
m.comparison <- function(data, export = NULL) {
  inertie <- ggplot(data) +
    geom_line(aes(m, Explained.inertia)) +
    geom_point(aes(m, Explained.inertia), color = "red") +
    labs(title = "a. Variation des données expliquées",
         y = "Inertie expliquée", x = "Paramètre m") + theme_tbe_gg()
  silhouette <- ggplot(data) +
    geom_line(aes(m, Silhouette.index)) +
    geom_point(aes(m, Silhouette.index), color = "red") +
    labs(title = "b. Consistance des groupes",
         y = "Critère de silhouette floue", x = "Paramètre m") + theme_tbe_gg()
  p <- ggarrange(inertie, silhouette)
  export.graphique(p, export); p
}

#' Optimisation du paramètre alpha (geocmeans) : silhouette et incohérence
alpha.comparison <- function(data, export = NULL) {
  silhouette <- ggplot(data) +
    geom_line(aes(alpha, Silhouette.index)) +
    geom_point(aes(alpha, Silhouette.index), color = "red") +
    labs(x = "Alpha", y = "Indice de silhouette") + theme_tbe_gg()
  incoherence <- ggplot(data) +
    geom_line(aes(alpha, spConsistency)) +
    geom_point(aes(alpha, spConsistency), color = "red") +
    labs(x = "Alpha", y = "Indice d'incohérence spatiale") + theme_tbe_gg()
  p <- ggarrange(silhouette, incoherence, ncol = 2, nrow = 1)
  export.graphique(p, export); p
}

#' Évolution temporelle du niveau moyen par groupe (geocmeans)
group.mean.comparison <- function(data, group.fields, export = NULL) {
  data_long <- pivot_longer(data, cols = all_of(group.fields),
                            names_to = "Cluster", values_to = "valeur")
  p <- ggplot(data_long, aes(year, valeur, color = Cluster, group = Cluster)) +
    geom_line() + geom_point() +
    labs(x = "Année", y = "Niveau TBE") +
    scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "yellow")) +
    theme_tbe_gg()
  export.graphique(p, export); p
}


# 5. GRAPHIQUES D'ÉVALUATION DES COVARIABLES --------------------------------

#' Classement des covariables par force d'association (barres horizontales)
plot.ranking <- function(ranking, outdir) {
  ranking$covariable <- factor(ranking$covariable, levels = rev(ranking$covariable))
  p <- ggplot(ranking, aes(x = covariable, y = force, fill = type)) +
    geom_col() + coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Force d'association avec l'intensité TBE",
         subtitle = "Spearman |rho| (continu) ou V de Cramér (catégoriel), 0 = nul, 1 = fort",
         x = NULL, y = "Force d'association [0-1]") + theme_tbe_gg()
  export.graphique(p, file.path(outdir, "classement_covariables.png"))
}

#' Heatmap de colinéarité entre covariables continues (corrplot)
plot.collinearity <- function(cor_mat, outdir) {
  if (!requireNamespace("corrplot", quietly = TRUE)) return(invisible())
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(outdir, "colinearite_covariables.png"), width = 700, height = 700)
  corrplot::corrplot(cor_mat, method = "color", type = "upper",
                     addCoef.col = "black", tl.col = "black",
                     title = "Colinéarité (Spearman) entre covariables continues",
                     mar = c(0, 0, 2, 0))
  dev.off()
}


# 6. GRAPHIQUES D'ÉVALUATION DES MODÈLES -----------------------------------

#' Matrice de confusion (heatmap : comptes + % normalisé par ligne)
save.confusion <- function(cm, model_name, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
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
         subtitle = "ligne = observé, colonne = prédit ; % par ligne (rappel en diagonale)",
         fill = "% ligne") + theme_tbe_gg(base_size = 11) +
    theme(legend.position = "right")
  export.graphique(p, file.path(outdir, paste0("confusion_", slug, ".png")),
                   width = 6.5, height = 5.5)
}

#' Comparaison du kappa pondéré entre modèles (CV spatiale vs test temporel)
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
         x = NULL, y = "Kappa pondéré (quadratique)", fill = NULL) + theme_tbe_gg()
  export.graphique(p, file.path(outdir, "kappa_comparaison.png"), width = 7.5, height = 5)
}

#' Rappel par classe d'intensité (détection des classes rares)
plot.class.recall <- function(class_tab, outdir, labels = ETATS_TBE) {
  df <- class_tab
  df$classe <- factor(df$classe, levels = labels)
  df$modele <- factor(df$modele, levels = unique(df$modele))
  p <- ggplot(df, aes(classe, rappel, fill = modele)) +
    geom_col(position = position_dodge(0.8), width = 0.7, na.rm = TRUE) +
    geom_text(aes(label = ifelse(is.na(rappel), "0", sprintf("%.2f", rappel))),
              position = position_dodge(0.8), vjust = -0.4, size = 2.8, na.rm = TRUE) +
    scale_fill_brewer(palette = "Set2") + ylim(0, 1) +
    labs(title = "Rappel par classe d'intensité (test temporel)",
         subtitle = "Part des cas réels de chaque état correctement prédits (1 = parfait)",
         x = NULL, y = "Rappel", fill = NULL) + theme_tbe_gg()
  export.graphique(p, file.path(outdir, "rappel_par_classe.png"), width = 8, height = 5)
}
