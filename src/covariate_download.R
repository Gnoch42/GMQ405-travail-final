# =============================================================================
# src/covariate_download.R
# -----------------------------------------------------------------------------
# TÉLÉCHARGEMENT CIBLÉ des données sources de covariables.
#
# Ne télécharge QUE ce qui intersecte la zone d'étude (study_zone), et SEULEMENT
# si ce n'est pas déjà présent localement (idempotent : relancer ne re-télécharge
# rien tant que la zone n'a pas changé). Les données écoforestières sont très
# volumineuses (~1 Go zippé par feuillet, ~3,6 Go décompressé) : on les découpe
# immédiatement à la zone d'étude et on supprime le brut (sauf keep_raw: true).
#
# Deux sources :
#   - Écoforestière MFFP : index de feuillets 250K -> ZIP GPKG par feuillet.
#   - Climat RSCQ        : index CSV -> fichiers NetCDF quotidiens par année.
#
# Usage :  Rscript src/covariate_download.R
#
# Dépendances : sf, yaml, utils (unzip/download.file), + src/features.R
# =============================================================================

library(sf)
library(yaml)
source("src/features.R")   # load.study.zone()

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Téléchargement ROBUSTE d'un gros fichier (curl : reprise + réessais)
#'
#' download.file() de base est peu fiable pour les fichiers volumineux (~1 Go)
#' du serveur MFFP : délai court, troncature silencieuse. On délègue à `curl`
#' avec reprise (-C -), réessais et suivi des redirections, puis on vérifie que
#' la taille finale correspond bien à l'en-tête du serveur.
robust.download <- function(url, dest) {
  status <- system2("curl", c("-fL", "-C", "-", "--retry", "4", "--retry-delay", "5",
                              "-o", shQuote(dest), shQuote(url)))
  if (status != 0 || !file.exists(dest))
    stop("Échec du téléchargement : ", url, " (code curl ", status, ")")
  # Vérification de complétude via l'en-tête Content-Length.
  hdr <- suppressWarnings(system2("curl", c("-sIL", shQuote(url)), stdout = TRUE))
  cl  <- grep("(?i)content-length", hdr, value = TRUE, perl = TRUE)
  if (length(cl) > 0) {
    expected <- as.numeric(sub("\\D+", "", tail(cl, 1)))
    if (!is.na(expected) && file.info(dest)$size < expected)
      stop("Fichier incomplet (", file.info(dest)$size, " < ", expected, " o) : ", url)
  }
  invisible(dest)
}


# -----------------------------------------------------------------------------
# 1. ÉCOFORESTIÈRE — feuillets intersectant la zone
# -----------------------------------------------------------------------------

#' Télécharger et découper les feuillets écoforestiers de la zone d'étude
#'
#' Pour chaque feuillet qui intersecte study_zone :
#'   - si la GPKG découpée existe déjà -> on saute (idempotent) ;
#'   - sinon : télécharge le ZIP (s'il manque), décompresse, découpe la couche
#'     des peuplements ET la table des essences à la zone, écrit une GPKG légère
#'     + un RDS d'essences, puis supprime le brut (sauf keep_raw).
#'
#' @return Vecteur des chemins des GPKG découpées produites.
download.ecoforest.tiles <- function(cfg) {
  eco <- cfg$covariate_sources$ecoforest
  target_crs <- cfg$project$target_crs %||% 32198
  zone <- load.study.zone(cfg$study_zone, target_crs)
  if (is.null(zone)) stop("study_zone requis pour le téléchargement écoforestier.")

  dir.create(eco$download_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(eco$clip_dir,     recursive = TRUE, showWarnings = FALSE)

  # Feuillets qui intersectent la zone (l'index est en WGS84).
  idx <- st_read(eco$index, quiet = TRUE)
  zone_idx <- st_transform(zone, st_crs(idx))
  suppressMessages({
    old_s2 <- sf_use_s2(FALSE)
    hit <- idx[st_intersects(idx, zone_idx, sparse = FALSE)[, 1], ]
    sf_use_s2(old_s2)
  })
  message("Feuillets écoforestiers intersectant la zone : ", nrow(hit),
          " (", paste(hit[[eco$tile_field]], collapse = ", "), ")")

  produced <- character(0)
  for (i in seq_len(nrow(hit))) {
    tile <- hit[[eco$tile_field]][i]
    url  <- hit[[eco$url_field]][i]
    clip_gpkg <- file.path(eco$clip_dir, paste0("ecoforest_", tile, ".gpkg"))

    # Idempotence : déjà découpé -> on saute complètement.
    if (file.exists(clip_gpkg) && !isTRUE(cfg$covariates_build$force)) {
      message("  [", tile, "] déjà découpé -> ignoré."); produced <- c(produced, clip_gpkg); next
    }

    zip_path <- file.path(eco$download_dir, basename(url))
    if (!file.exists(zip_path)) {
      message("  [", tile, "] téléchargement du ZIP (~1 Go)...")
      robust.download(url, zip_path)
    }

    # Sélection du bon .gpkg dans le ZIP. Certains feuillets contiennent DEUX
    # variantes dans un sous-dossier : "_C" (compilée, avec attributs dérivés) et
    # "_NC" (non compilée). On privilégie la variante compilée ; sinon l'unique
    # .gpkg présent.
    entries <- utils::unzip(zip_path, list = TRUE)$Name
    gpkgs   <- entries[grepl("\\.gpkg$", entries, ignore.case = TRUE)]
    compiled <- gpkgs[grepl("_C\\.gpkg$", gpkgs, ignore.case = TRUE)]
    gpkg_name <- if (length(compiled) > 0) compiled[1] else gpkgs[1]
    raw_gpkg <- file.path(eco$download_dir, gpkg_name)
    if (!file.exists(raw_gpkg)) {
      message("  [", tile, "] décompression (", basename(gpkg_name), ")...")
      utils::unzip(zip_path, files = gpkg_name, exdir = eco$download_dir)
    }

    clip.ecoforest.tile(raw_gpkg, zone, eco, clip_gpkg, tile)
    produced <- c(produced, clip_gpkg)

    # Nettoyage du brut (volumineux) sauf demande contraire.
    if (!isTRUE(eco$keep_raw)) {
      unlink(c(zip_path, raw_gpkg))
      message("  [", tile, "] brut supprimé (keep_raw = false).")
    }
  }
  produced
}


#' Découper une GPKG écoforestière brute à la zone d'étude
#'
#' Écrit deux fichiers légers : la couche des peuplements clippée (GPKG) et la
#' table des essences correspondantes (RDS, jointe par geocode).
clip.ecoforest.tile <- function(raw_gpkg, zone, eco, clip_gpkg, tile) {
  wkt <- st_as_text(st_geometry(st_transform(zone, 32198)))

  # Peuplements : lecture filtrée spatialement (rapide via l'index de la GPKG).
  suppressMessages(sf_use_s2(FALSE))
  pee <- st_read(raw_gpkg, layer = eco$stand_layer, wkt_filter = wkt, quiet = TRUE)
  keep_cols <- intersect(c("geocode", "gr_ess", "cl_age", "type_couv", "superficie"),
                         names(pee))
  pee <- pee[, keep_cols]
  st_write(pee, clip_gpkg, layer = eco$stand_layer, delete_dsn = TRUE, quiet = TRUE)

  # Table des essences (non spatiale) : on ne garde que les geocodes clippés.
  ess <- st_read(raw_gpkg,
                 query = sprintf("SELECT geocode, etage, essence, st_ess_pc FROM %s",
                                 eco$essence_table), quiet = TRUE)
  ess <- ess[ess$geocode %in% pee$geocode, ]
  saveRDS(ess, sub("\\.gpkg$", "_essence.rds", clip_gpkg))

  message("  [", tile, "] découpé : ", nrow(pee), " peuplements, ",
          nrow(ess), " lignes d'essences.")
}


# -----------------------------------------------------------------------------
# 2. CLIMAT RSCQ — fichiers NetCDF par année
# -----------------------------------------------------------------------------

#' Télécharger les NetCDF climatiques pour un paramètre et des années données
#'
#' @param years Années nécessaires (inclure t-1 si l'hiver déborde de décembre).
#' @return Vecteur des chemins des fichiers .nc disponibles localement.
download.climate.files <- function(cfg, parameter, years) {
  clim <- cfg$covariate_sources$climate
  dir.create(clim$download_dir, recursive = TRUE, showWarnings = FALSE)

  links <- utils::read.csv(clim$links_csv, stringsAsFactors = FALSE)
  sel <- links[links[[clim$param_field]] == parameter &
               links[[clim$year_field]] %in% years, ]
  if (nrow(sel) == 0) {
    warning("Aucun fichier climatique trouvé pour ", parameter, " / années demandées.")
    return(character(0))
  }

  paths <- character(0)
  for (i in seq_len(nrow(sel))) {
    fname <- sel[[clim$file_field]][i]
    dest  <- file.path(clim$download_dir, fname)
    if (!file.exists(dest)) {
      message("  climat : téléchargement ", fname, " (~",
              round(sel$Taille_Ko[i] / 1024), " Mo)...")
      robust.download(sel[[clim$url_field]][i], dest)
    } else {
      message("  climat : ", fname, " déjà présent -> ignoré.")
    }
    paths <- c(paths, dest)
  }
  paths
}


# -----------------------------------------------------------------------------
# ORCHESTRATION
# -----------------------------------------------------------------------------

#' Télécharger toutes les sources nécessaires selon la config
main.download <- function(config_path = "config.yaml") {
  cfg <- read_yaml(config_path)

  message("== Téléchargement écoforestier ==")
  download.ecoforest.tiles(cfg)

  # Années climatiques nécessaires : celles de la cible + l'année précédente
  # (l'hiver de l'année t inclut décembre de t-1).
  message("\n== Téléchargement climatique (TMIN) ==")
  tgt <- st_read(cfg$target$source, quiet = TRUE)
  yrs <- sort(unique(tgt[[cfg$target$year_field]]))
  years_needed <- sort(unique(c(yrs, yrs - 1)))
  param <- cfg$covariates_build$items$tmin_hiver$parameter %||% "TMIN"
  download.climate.files(cfg, param, years_needed)

  message("\nTéléchargements terminés.")
}

if (sys.nframe() == 0) {
  main.download("config.yaml")
}
