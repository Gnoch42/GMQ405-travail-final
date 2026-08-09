library(sf)
source("src/utils.R")   # robust.download(), %||%

extract.tbe.data <- function(cfg) {
  tbe <- st_read(cfg$extraction$TBE$inpath)
  sda <- st_read(cfg$extraction$SDA$inpath, layer = cfg$extraction$SDA$layer)
  
  region.of.interest <- subset(sda, st_drop_geometry(sda)[,cfg$extraction$SDA$field] %in% cfg$extraction$SDA$regions) |>
    st_union()
  if (st_crs(region.of.interest) != st_crs(tbe)) {
    region.of.interest <- st_transform(region.of.interest, st_crs(tbe))
  }
  tbe.subset <- subset(tbe, st_intersects(tbe, region.of.interest, sparse = FALSE))
  
  st_write(tbe.subset, cfg$extraction$TBE$outpath, append = FALSE) 
}


extract.ieqm.data <- function(cfg) {
  ieqm.list <- list.files(cfg$extraction$IEQM$inpath, full.names = TRUE)
  sda <- st_read(cfg$extraction$SDA$inpath, layer = cfg$extraction$SDA$layer)
  
  region.of.interest <- subset(sda, st_drop_geometry(sda)[,cfg$extraction$SDA$field] %in% cfg$extraction$SDA$regions) |>
    st_union()
  output.data <- NULL
  for(file in ieqm.list){
    ieqm <- st_read(file, layer = 'pee_ori')
    if (st_crs(region.of.interest) != st_crs(ieqm)) {
      region.of.interest <- st_transform(region.of.interest, st_crs(ieqm))
    }
    ieqm.subset <- subset(ieqm, st_intersects(region.of.interest, ieqm, sparse = FALSE))
    if(is.null(output.data)){
      output.data <- ieqm.subset
    } else {
      output.data <- rbind(output.data, ieqm.subset)
    }
  }
  output.data <- subset(output.data, !duplicated(output.data))
  st_write(output.data, cfg$extraction$IEQM$outpath, append = FALSE)
}


#' Emprise de la zone d'étude (union des MRC), dans un CRS donné
.study.region <- function(cfg, crs_target) {
  sda <- st_read(cfg$study_zone$source, layer = cfg$study_zone$layer, quiet = TRUE)
  region <- subset(sda, st_drop_geometry(sda)[, cfg$study_zone$field] %in% cfg$study_zone$regions)
  st_transform(st_union(st_zm(region)), crs_target)
}

#' Télécharger les blocs GRHQ intersectant la zone (idempotent)
#'
#' Utilise l'INDEX des blocs (champ URL) : ne télécharge que les blocs qui
#' touchent la zone d'étude, et seulement si le FGDB correspondant n'est pas déjà
#' présent. Chaque bloc est décompressé dans son propre sous-dossier.
#'
#' @return Vecteur des chemins des géodatabases (.gdb) disponibles localement.
download.grhq.blocs <- function(cfg) {
  g <- cfg$extraction$GRHQ
  dir.create(g$download_dir, recursive = TRUE, showWarnings = FALSE)

  idx  <- st_read(g$index, quiet = TRUE)
  zone <- .study.region(cfg, st_crs(idx))
  old  <- sf::sf_use_s2(FALSE)
  hit  <- idx[st_intersects(idx, zone, sparse = FALSE)[, 1], ]
  sf::sf_use_s2(old)
  message("Blocs GRHQ intersectant la zone : ", nrow(hit), " (",
          paste(hit[[g$bloc_field]], collapse = ", "), ")")

  fgdbs <- character(0)
  for (i in seq_len(nrow(hit))) {
    bloc <- hit[[g$bloc_field]][i]
    url  <- hit[[g$url_field]][i]
    bloc_dir <- file.path(g$download_dir, bloc)
    gdb <- list.files(bloc_dir, pattern = "\\.gdb$", full.names = TRUE, include.dirs = TRUE)

    if (length(gdb) == 0 || isTRUE(g$force_extraction)) {
      dir.create(bloc_dir, recursive = TRUE, showWarnings = FALSE)
      zip <- file.path(bloc_dir, basename(url))
      if (!file.exists(zip)) {
        message("  [", bloc, "] téléchargement du ZIP...")
        robust.download(url, zip)
      }
      message("  [", bloc, "] décompression...")
      utils::unzip(zip, exdir = bloc_dir)
      gdb <- list.files(bloc_dir, pattern = "\\.gdb$", full.names = TRUE, include.dirs = TRUE)
    } else {
      message("  [", bloc, "] déjà présent -> ignoré.")
    }
    fgdbs <- c(fgdbs, gdb)
  }
  fgdbs
}

#' Extraire le fleuve Saint-Laurent depuis les blocs GRHQ de la zone
#'
#' Télécharge les blocs nécessaires via l'index (si absents), lit la couche du
#' réseau hydrographique surfacique, ne garde que le fleuve, puis fusionne.
extract.St_Laurent <- function(cfg) {
  g <- cfg$extraction$GRHQ
  fgdbs <- download.grhq.blocs(cfg)

  layer      <- g$layer          %||% "RH_S"
  topo_field <- g$toponyme_field %||% "TOPONYME"
  topo       <- g$toponyme       %||% "Fleuve Saint-Laurent"

  output.data <- NULL
  for (fgdb in fgdbs) {
    # Requête SQL (filtre à la lecture) avec repli sur lecture + filtre en R.
    q <- sprintf("SELECT * FROM %s WHERE %s = '%s'", layer, topo_field, topo)
    grhq <- tryCatch(st_read(fgdb, query = q, quiet = TRUE),
                     error = function(e) {
                       r <- st_read(fgdb, layer = layer, quiet = TRUE)
                       r[st_drop_geometry(r)[[topo_field]] == topo, ]
                     })
    if (nrow(grhq) == 0) next
    output.data <- if (is.null(output.data)) grhq else rbind(output.data, grhq)
  }

  output.data <- subset(output.data, !duplicated(output.data)) |>
    st_union() |>
    st_as_sf()
  dir.create(dirname(g$outpath), recursive = TRUE, showWarnings = FALSE)
  st_write(output.data, g$outpath, append = FALSE, quiet = TRUE)
}
