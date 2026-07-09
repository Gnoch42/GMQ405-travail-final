# Fonction pour créer une grille hexagonale par-dessus la zone d'étude

#' Créer une grille hexagonale couvrant la zone d'étude
#'
#' @param cfg Liste de configuration lue depuis config.yaml
#'
#' @return Un objet sf contenant la grille hexagonale (invisible),
#'   également écrit sur disque au chemin cfg$hex_grid$outpath
#'
#' @details
#'   La zone d'étude est dérivée des données IEQM déjà extraites
#'   (cfg$extraction$IEQM$outpath). La taille des cellules hexagonales
#'   est lue depuis cfg$hex_grid$cellsize (en unités du CRS, typiquement
#'   des mètres si le CRS est projeté). Seuls les hexagones qui
#'   intersectent réellement la zone d'étude sont conservés.
create.hex.grid <- function(cfg) {

  # --- 1. Chargement de la zone d'étude ---
  # On utilise les données IEQM extraites comme emprise de la zone d'étude.
  # st_union() fusionne tous les polygones/points en une seule géométrie
  # qui représente l'enveloppe globale de la zone d'intérêt.
  ieqm.data <- st_read(cfg$extraction$IEQM$outpath, quiet = TRUE)
  zone <- st_union(ieqm.data)

  # --- 2. Récupération de la taille de cellule depuis la config ---
  cellsize <- cfg$hex_grid$cellsize

  # --- 3. Création de la grille hexagonale ---
  # st_make_grid() génère une grille couvrant la bounding box de `zone`.
  # square = FALSE produit des hexagones au lieu de carrés.
  # flat_topped = FALSE : hexagones pointus vers le haut (orientation par défaut).
  hex.grid <- st_make_grid(
    zone,
    cellsize  = cellsize,
    square    = FALSE,
    flat_topped = FALSE
  )

  # Conversion en objet sf et ajout d'un identifiant unique par hexagone
  hex.grid <- st_sf(id = seq_along(hex.grid), geometry = hex.grid)

  # --- 4. Découpage à la zone d'étude ---
  # On ne conserve que les hexagones qui intersectent réellement la zone
  # d'étude, afin d'éviter une grille trop large couvrant des zones vides.
  intersects.zone <- st_intersects(hex.grid, zone, sparse = FALSE)[, 1]
  hex.grid <- hex.grid[intersects.zone, ]

  message("Grille hexagonale créée : ", nrow(hex.grid), " hexagones (cellsize = ", cellsize, ")")

  # --- 5. Écriture du résultat sur disque ---
  st_write(hex.grid, cfg$hex_grid$outpath, append = FALSE, quiet = TRUE)

  invisible(hex.grid)
}
