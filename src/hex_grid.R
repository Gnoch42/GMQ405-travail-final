# =============================================================================
# src/hex_grid.R
# -----------------------------------------------------------------------------
# Grille hexagonale + agrégation spatiale générique.
#
# Ce fichier fournit les DEUX briques de base réutilisées partout dans le
# projet (modélisation ET évaluation des covariables) :
#
#   1. create.hex.grid()  : fabrique la grille hexagonale régulière qui sert
#                           d'unité spatiale commune.
#   2. aggregate.to.hex() : ramène N'IMPORTE QUELLE source de données
#                           (polygones OU raster, de forme quelconque) sur
#                           cette grille, en assignant une valeur à chaque
#                           hexagone selon ce qu'il contient.
#
# Pourquoi une grille hexagonale ?
#   Nos données sources (intensité TBE, covariables) arrivent sous des formes
#   très variées (polygones irréguliers, rasters de résolutions différentes).
#   Pour comparer et modéliser, il faut une unité spatiale UNIQUE et régulière.
#   L'hexagone est un bon compromis : contrairement au carré, tous les voisins
#   d'un hexagone sont à la même distance de son centre, ce qui reflète mieux
#   les phénomènes de propagation (comme la progression d'une infestation).
#
# Dépendances : sf (vecteur), terra (raster).
# =============================================================================

library(sf)
library(terra)


## HELPERS D'AGRÉGATION PONDÉRÉE ------------------------------------------
# Ces deux fonctions résument un ensemble de valeurs en UNE seule valeur pour
# un hexagone, en tenant compte du POIDS de chaque valeur (surface d'intersection
# pour les polygones, fraction de cellule couverte pour les rasters).
#
# Pourquoi pondérer ? Si un hexagone est recouvert à 90 % par un polygone
# "forêt dense" et à 10 % par un polygone "clairière", un simple décompte
# (1 vs 1) donnerait une égalité trompeuse. La pondération par surface reflète
# ce qui occupe RÉELLEMENT l'hexagone. Ce choix change les résultats : il est
# donc explicite et documenté ici.

#' Médiane pondérée
#'
#' Utilisée pour les variables QUANTITATIVES CONTINUES (ex. température,
#' altitude). La médiane (plutôt que la moyenne) est robuste aux valeurs
#' extrêmes : une cellule aberrante ne tire pas tout l'hexagone vers le haut.
#'
#' @param x Vecteur numérique de valeurs.
#' @param w Vecteur de poids (même longueur que x). Par défaut, poids égaux.
#' @return La valeur v telle que la moitié du poids total se trouve de part
#'   et d'autre de v.
weighted.median <- function(x, w = rep(1, length(x))) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  x <- x[keep]; w <- w[keep]
  if (length(x) == 0) return(NA_real_)
  # Cas dégénérés : une seule valeur (ex. hexagone couvrant une seule cellule
  # raster) ou valeurs toutes identiques -> pas d'interpolation possible/utile.
  if (length(x) == 1 || length(unique(x)) == 1) return(x[1])

  ord <- order(x)                       # tri des valeurs par ordre croissant
  x <- x[ord]; w <- w[ord]
  cum <- cumsum(w) - 0.5 * w            # poids cumulé "centré" sur chaque valeur
  # On interpole pour trouver la valeur au point où 50 % du poids est atteint.
  approx(cum, x, xout = 0.5 * sum(w), rule = 2, ties = "ordered")$y
}


#' Mode pondéré (catégorie dominante)
#'
#' Utilisé pour les variables QUALITATIVES (catégorielles ou ordinales, ex.
#' essence dominante, classe d'intensité TBE). Renvoie la catégorie qui
#' occupe la plus grande surface / le plus grand nombre de cellules dans
#' l'hexagone.
#'
#' @param x Vecteur de valeurs (facteur, caractère ou entier).
#' @param w Vecteur de poids (surface d'intersection ou couverture cellule).
#' @return La catégorie de poids total maximal (même type que x en entrée).
weighted.mode <- function(x, w = rep(1, length(x))) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  x <- x[keep]; w <- w[keep]
  if (length(x) == 0) return(NA)

  # Somme des poids par catégorie, puis on garde la catégorie au poids max.
  totals <- tapply(w, x, sum)
  winner <- names(totals)[which.max(totals)]

  # On restitue le type d'origine (un facteur reste un facteur, etc.).
  if (is.factor(x))      factor(winner, levels = levels(x))
  else if (is.numeric(x)) as(winner, class(x))
  else                    winner
}


## 1. CRÉATION DE LA GRILLE HEXAGONALE --------------------------------------

#' Créer une grille hexagonale couvrant une zone d'étude
#'
#' @param zone       Objet sf/sfc délimitant la zone d'étude (ex. l'emprise des
#'                   données TBE, ou un contour administratif). La grille ne
#'                   conserve que les hexagones intersectant cette zone.
#' @param cellsize   Taille de la cellule, exprimée dans les UNITÉS DU CRS cible
#'                   (voir target_crs). Avec un CRS projeté en mètres, cellsize
#'                   correspond à la distance entre deux côtés opposés de
#'                   l'hexagone (ex. 5000 = hexagones d'environ 5 km).
#' @param target_crs (Optionnel) CRS projeté cible (code EPSG ou objet crs).
#'                   IMPORTANT : la grille DOIT être construite dans un CRS
#'                   projeté (en mètres) pour que `cellsize` ait un sens métrique.
#'                   Les données TBE sont en NAD83 géographique (EPSG:4269, en
#'                   degrés) ; on reprojette donc vers, par défaut, le Lambert
#'                   québécois (EPSG:32198, en mètres). Si NULL, on garde le CRS
#'                   de `zone` tel quel (à n'utiliser que s'il est déjà projeté).
#'
#' @return Un objet sf : une ligne par hexagone, colonnes `hex_id` (identifiant
#'   unique), `hex_area` (surface de l'hexagone) et `geometry`.
create.hex.grid <- function(zone, cellsize, target_crs = 32198) {

  # --- 1. Reprojection vers un CRS projeté (mètres) ---
  # Sans cela, cellsize serait interprété en degrés, ce qui n'a aucun sens
  # comme distance et déformerait fortement les hexagones selon la latitude.
  if (!is.null(target_crs)) {
    zone <- st_transform(zone, target_crs)
  }

  # --- 2. Fusion de la zone en une seule emprise ---
  # st_union() combine tous les polygones en une géométrie unique servant de
  # gabarit (on ne veut pas une grille par polygone).
  zone <- st_union(zone)

  # --- 3. Génération de la grille hexagonale ---
  # st_make_grid() pave la bounding box de la zone.
  #   square = FALSE     -> hexagones (au lieu de carrés)
  #   flat_topped = FALSE -> hexagones "pointe en haut" (orientation par défaut)
  hex <- st_make_grid(zone, cellsize = cellsize, square = FALSE, flat_topped = FALSE)
  hex <- st_sf(hex_id = seq_along(hex), geometry = hex)

  # --- 4. Découpage à la zone d'étude ---
  # On ne garde que les hexagones qui touchent réellement la zone, pour éviter
  # une grille rectangulaire débordant sur des zones sans données.
  touches <- st_intersects(hex, zone, sparse = FALSE)[, 1]
  hex <- hex[touches, ]

  # --- 5. Ré-indexation + surface de référence ---
  # hex_area sert ensuite à calculer le taux de couverture de chaque hexagone
  # (pour gérer les hexagones de bordure partiellement couverts).
  hex$hex_id   <- seq_len(nrow(hex))
  hex$hex_area <- as.numeric(st_area(hex))

  message("Grille hexagonale : ", nrow(hex), " hexagones (cellsize = ", cellsize,
          ", CRS = ", st_crs(hex)$epsg, ")")
  hex
}


## 2. AGRÉGATION D'UNE SOURCE (VECTEUR OU RASTER) SUR LA GRILLE -------------

#' Assigner une valeur à chaque hexagone à partir d'une source quelconque
#'
#' Fonction PIVOT du projet : elle prend une source de données de forme
#' arbitraire (polygones OU raster) et calcule, pour chaque hexagone, une valeur
#' résumée de ce que la source contient à cet endroit.
#'
#' Stratégie selon le type de source :
#'   - VECTEUR (polygones) : intersection spatiale hexagone × polygones, puis
#'     agrégation pondérée par la SURFACE d'intersection.
#'   - RASTER : extraction zonale (quelles cellules tombent dans l'hexagone),
#'     puis agrégation pondérée par la FRACTION de chaque cellule couverte.
#'
#' Stratégie selon le type de variable :
#'   - continue  -> médiane pondérée  (weighted.median)
#'   - qualitative (catégorielle/ordinale) -> mode pondéré (weighted.mode)
#'
#' Gestion des bordures : un hexagone dont la source ne couvre qu'une petite
#' fraction (< min_coverage) reçoit NA, pour éviter des valeurs peu fiables
#' basées sur un échantillon minuscule.
#'
#' @param hex          Grille hexagonale (sortie de create.hex.grid()).
#' @param source       Source de données : objet sf/sfc (vecteur) OU SpatRaster
#'                     (terra). Le type est détecté automatiquement.
#' @param value_field  Pour un VECTEUR : nom de la colonne contenant la valeur.
#'                     Ignoré pour un raster (on utilise la bande `band`).
#' @param value_type   "continuous" (médiane) ou "categorical"/"ordinal" (mode).
#' @param out_name     Nom de la colonne de sortie ajoutée à la grille.
#' @param min_coverage Fraction minimale de l'hexagone devant être couverte par
#'                     la source pour produire une valeur (défaut 0.25 = 25 %).
#'                     En dessous, la valeur est NA.
#' @param band         Pour un RASTER : indice de la bande à extraire (défaut 1).
#' @param add_coverage Si TRUE, ajoute une colonne `<out_name>_cov` avec le taux
#'                     de couverture (utile pour diagnostiquer les bordures).
#'
#' @return La grille `hex` enrichie d'une colonne `out_name` (et éventuellement
#'   de la couverture).
aggregate.to.hex <- function(hex, source, value_field = NULL,
                             value_type = c("continuous", "categorical", "ordinal"),
                             out_name = "value", min_coverage = 0.25,
                             band = 1, add_coverage = FALSE) {

  value_type <- match.arg(value_type)
  is_continuous <- value_type == "continuous"

  # Vecteur de résultats (un par hexagone) + couverture, initialisés à NA/0.
  n <- nrow(hex)
  result   <- vector("list", n)   # liste pour préserver les types (facteur, etc.)
  coverage <- numeric(n)

  # ===========================================================================
  # CAS 1 — SOURCE RASTER (SpatRaster de terra)
  # ===========================================================================
  if (inherits(source, "SpatRaster")) {

    # Reprojection du raster vers le CRS de la grille si nécessaire.
    if (!crs.matches(source, hex)) {
      # 'near' (plus proche voisin) pour les catégories afin de ne pas inventer
      # de valeurs intermédiaires ; 'bilinear' pour le continu.
      source <- project(source, paste0("EPSG:", st_crs(hex)$epsg),
                        method = if (is_continuous) "bilinear" else "near")
    }

    hex_v    <- vect(hex)                              # sf -> SpatVector (terra)
    cell_area <- prod(res(source))                     # surface d'une cellule
    band_lyr <- source[[band]]

    # terra::extract(..., weights=TRUE) renvoie, pour chaque hexagone, la liste
    # des cellules intersectées avec `weight` = fraction de la cellule couverte.
    ex <- terra::extract(band_lyr, hex_v, weights = TRUE, ID = TRUE)
    colnames(ex) <- c("ID", "val", "weight")

    for (i in seq_len(n)) {
      rows <- ex[ex$ID == i, , drop = FALSE]
      if (nrow(rows) == 0) { result[[i]] <- NA; next }

      # Couverture = surface de cellules couvertes / surface de l'hexagone.
      coverage[i] <- sum(rows$weight, na.rm = TRUE) * cell_area / hex$hex_area[i]
      if (coverage[i] < min_coverage) { result[[i]] <- NA; next }

      result[[i]] <- if (is_continuous) {
        weighted.median(rows$val, rows$weight)
      } else {
        weighted.mode(rows$val, rows$weight)
      }
    }

  # ===========================================================================
  # CAS 2 — SOURCE VECTEUR (polygones sf/sfc)
  # ===========================================================================
  } else if (inherits(source, c("sf", "sfc"))) {

    if (is.null(value_field) && inherits(source, "sf"))
      stop("value_field doit être fourni pour une source vecteur.")

    source <- st_transform(source, st_crs(hex))       # aligner les CRS

    # Intersection : découpe chaque polygone source par chaque hexagone.
    # st_intersection conserve les attributs de la source ET l'identité de
    # l'hexagone (hex_id). Chaque ligne = un morceau (hexagone × polygone).
    # warnings d'attributs supposés constants : sans conséquence ici.
    pieces <- suppressWarnings(st_intersection(hex[, c("hex_id", "hex_area")], source))
    if (nrow(pieces) > 0) {
      pieces$piece_area <- as.numeric(st_area(pieces))
      vals   <- st_drop_geometry(pieces)[[value_field]]
      hid    <- pieces$hex_id
      parea  <- pieces$piece_area
      harea  <- pieces$hex_area

      # Couverture de chaque hexagone = surface intersectée / surface hexagone.
      cov_by_hex <- tapply(parea, hid, sum) / tapply(harea, hid, function(a) a[1])

      for (i in seq_len(n)) {
        id <- hex$hex_id[i]
        sel <- hid == id
        if (!any(sel)) { result[[i]] <- NA; next }

        coverage[i] <- as.numeric(cov_by_hex[as.character(id)])
        if (is.na(coverage[i]) || coverage[i] < min_coverage) { result[[i]] <- NA; next }

        result[[i]] <- if (is_continuous) {
          weighted.median(vals[sel], parea[sel])
        } else {
          weighted.mode(vals[sel], parea[sel])
        }
      }
    }

  } else {
    stop("Type de source non reconnu : fournir un objet sf/sfc ou un SpatRaster.")
  }

  # --- Assemblage de la colonne de sortie ---
  # unlist() reconstitue un vecteur atomique en préservant le type dominant.
  out_vec <- if (is_continuous) {
    vapply(result, function(v) if (is.null(v)) NA_real_ else as.numeric(v), numeric(1))
  } else {
    # Pour le qualitatif on garde des libellés caractères (mode).
    vapply(result, function(v) if (length(v) == 0 || is.null(v)) NA_character_ else as.character(v), character(1))
  }

  hex[[out_name]] <- out_vec
  if (add_coverage) hex[[paste0(out_name, "_cov")]] <- coverage
  hex
}


## Petit utilitaire : comparaison de CRS entre un raster terra et un objet sf. ----
crs.matches <- function(rast, hex) {
  isTRUE(try(terra::same.crs(rast, paste0("EPSG:", st_crs(hex)$epsg)), silent = TRUE))
}
