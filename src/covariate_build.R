# =============================================================================
# src/covariate_build.R
# -----------------------------------------------------------------------------
# CONSTRUCTION des covariables réelles et assemblage du PANEL hex×année.
#
# Produit un fichier de CACHE unique (data/covariates/panel_hex.rds + .gpkg) qui
# contient, pour chaque hexagone et chaque année : la cible TBE et les 5
# covariables. Ce cache est ensuite lu SANS RECALCUL par :
#   - l'outil d'évaluation des covariables (src/covariate_evaluation.R),
#   - le pipeline de modélisation (mod_spa_temp.R).
#
# Les 5 covariables (définies dans config `covariates_build$items`) :
#   prop_hote       (statique)  proportion surfacique sapin+épinettes / forêt
#   age_peuplement  (statique)  âge médian du peuplement (étage dominant)
#   tmin_hiver      (annuelle)  min de TMIN quotidien, Déc(t-1)+Jan-Fév(t)
#   hist_infest     (annuelle)  somme des sévérités TBE sur t-5..t-1
#   dist_foyer      (annuelle)  distance au foyer actif le plus proche à t-1
#
# "statique" = une valeur par hexagone (valable toutes années) ; "annuelle" =
# une valeur par hexagone ET par année.
#
# Usage :  Rscript src/covariate_build.R   (après covariate_download.R)
#
# Dépendances : sf, terra, yaml, + src/hex_grid.R, src/features.R
# =============================================================================

library(sf)
library(terra)
library(yaml)
source("src/hex_grid.R")
source("src/features.R")

`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# OUTILS ÉCOFORESTIERS
# -----------------------------------------------------------------------------

#' Charger les peuplements + essences découpés (tous feuillets de la zone)
load.ecoforest <- function(cfg) {
  eco <- cfg$covariate_sources$ecoforest
  gpkgs <- list.files(eco$clip_dir, pattern = "^ecoforest_.*\\.gpkg$", full.names = TRUE)
  if (length(gpkgs) == 0)
    stop("Aucune GPKG écoforestière découpée dans ", eco$clip_dir,
         " : lancez d'abord src/covariate_download.R")

  pee_list <- lapply(gpkgs, function(g) st_read(g, layer = eco$stand_layer, quiet = TRUE))
  ess_list <- lapply(sub("\\.gpkg$", "_essence.rds", gpkgs), readRDS)
  pee <- do.call(rbind, pee_list)
  ess <- do.call(rbind, ess_list)
  list(pee = st_transform(pee, cfg$project$target_crs %||% 32198), essence = ess)
}


#' Convertir un code de classe d'âge MFFP en âge numérique (années)
#'
#' Le champ `cl_age` est un CODE, pas un âge : classes équiennes ("10","30",...,
#' "120"), peuplements à deux étages ("3050" = 30 sur 50), et classes
#' INÉQUIENNES ("JIN"/"JIR" jeune, "VIN"/"VIR" vieux). On retient l'âge de
#' l'ÉTAGE DOMINANT (premier code lu), convention documentée et ajustable.
#'
#' @param codes Vecteur de codes cl_age.
#' @param age_ji Âge attribué aux classes jeunes inéquiennes (JI*).
#' @param age_vi Âge attribué aux classes vieilles inéquiennes (VI*).
#' @return Vecteur numérique d'âges (NA si non reconnu / null).
parse.age.class <- function(codes, age_ji = 30, age_vi = 90) {
  sapply(codes, function(cd) {
    if (is.na(cd) || cd == "") return(NA_real_)
    # Jeton de tête : 120 (3 car.), puis JI/VI (inéquienne), puis 2 chiffres.
    if (startsWith(cd, "120")) return(120)
    lead2 <- substr(cd, 1, 2)
    if (lead2 %in% c("JI")) return(age_ji)
    if (lead2 %in% c("VI")) return(age_vi)
    if (lead2 %in% c("10","30","50","70","90")) return(as.numeric(lead2))
    NA_real_
  }, USE.NAMES = FALSE)
}


# -----------------------------------------------------------------------------
# CONSTRUCTEURS DE COVARIABLES
# -----------------------------------------------------------------------------
# Chaque constructeur renvoie un data.frame :
#   - statique : colonnes (hex_id, <nom>)
#   - annuelle : colonnes (hex_id, year, <nom>)

#' prop_hote : proportion surfacique d'essences hôtes sur la forêt de l'hexagone
#'
#' Pour chaque peuplement, `host_frac` = part de surface terrière en essences
#' hôtes (somme des st_ess_pc hôtes / 100). Par hexagone :
#'   prop_hote = Σ(aire∩ × host_frac) / Σ(aire∩)   sur les peuplements FORESTIERS.
#' C'est la fraction de la forêt de l'hexagone occupée par sapin+épinettes.
build.prop.hote <- function(hex, eco_data, item, name) {
  pee <- eco_data$pee
  ess <- eco_data$essence

  # % hôte par geocode (somme des essences hôtes).
  host <- ess[ess$essence %in% item$host_essences, ]
  host_pc <- tapply(host$st_ess_pc, host$geocode, sum, na.rm = TRUE)
  pee$host_frac <- pmin(1, (host_pc[pee$geocode] %||% 0) / 100)
  pee$host_frac[is.na(pee$host_frac)] <- 0
  pee$is_forest <- !is.na(pee$type_couv) & pee$type_couv %in% item$forest_cover

  # Intersection peuplements × hexagones -> aires pondérées.
  pee_f <- pee[pee$is_forest, c("host_frac")]
  suppressWarnings({
    old <- sf_use_s2(FALSE)
    pieces <- st_intersection(hex[, "hex_id"], pee_f)
    sf_use_s2(old)
  })
  pieces$a <- as.numeric(st_area(pieces))
  d <- st_drop_geometry(pieces)
  num <- tapply(d$a * d$host_frac, d$hex_id, sum, na.rm = TRUE)
  den <- tapply(d$a,               d$hex_id, sum, na.rm = TRUE)
  val <- num / den

  out <- data.frame(hex_id = hex$hex_id)
  out[[name]] <- as.numeric(val[as.character(hex$hex_id)])
  out
}


#' age_peuplement : âge médian (surface-pondéré) de l'étage dominant
build.age <- function(hex, eco_data, item, name, min_coverage) {
  pee <- eco_data$pee
  pee$age_num <- parse.age.class(pee$cl_age,
                                 item$age_jeune_inequienne %||% 30,
                                 item$age_vieux_inequienne %||% 90)
  pee <- pee[!is.na(pee$age_num), "age_num"]
  h <- aggregate.to.hex(hex, pee, value_field = "age_num", value_type = "continuous",
                        out_name = name, min_coverage = min_coverage)
  st_drop_geometry(h)[, c("hex_id", name)]
}


#' tmin_hiver : min de TMIN quotidien sur Déc(t-1)+Jan-Fév(t), par année
#'
#' Lit les NetCDF (déjà en EPSG:32198, indexés par date via time()), extrait les
#' jours d'hiver, prend le minimum par cellule, puis agrège par la médiane sur
#' chaque hexagone. Renvoie une valeur par hexagone ET par année.
build.winter.tmin <- function(hex, cfg, item, name, years, min_coverage) {
  clim_dir <- cfg$covariate_sources$climate$download_dir
  param <- item$parameter %||% "TMIN"
  prev_m <- item$winter_prev_months %||% c(12)
  curr_m <- item$winter_curr_months %||% c(1, 2)

  nc_path <- function(y) file.path(clim_dir, sprintf("%s_%d.nc", param, y))

  rows <- list()
  for (t in years) {
    layers <- list()
    # Décembre de t-1.
    if (file.exists(nc_path(t - 1))) {
      r <- rast(nc_path(t - 1))
      m <- as.integer(format(terra::time(r), "%m"))
      if (any(m %in% prev_m)) layers[["prev"]] <- r[[which(m %in% prev_m)]]
    }
    # Jan-fév de t.
    if (file.exists(nc_path(t))) {
      r <- rast(nc_path(t))
      m <- as.integer(format(terra::time(r), "%m"))
      if (any(m %in% curr_m)) layers[["curr"]] <- r[[which(m %in% curr_m)]]
    }
    if (length(layers) == 0) {   # fichiers manquants -> NA
      df <- data.frame(hex_id = hex$hex_id, year = t); df[[name]] <- NA_real_
      rows[[as.character(t)]] <- df; next
    }
    winter <- do.call(c, unname(layers))
    wmin <- app(winter, "min", na.rm = TRUE)              # min hivernal par cellule
    h <- aggregate.to.hex(hex, wmin, value_type = "continuous",
                          out_name = name, min_coverage = min_coverage)
    df <- st_drop_geometry(h)[, c("hex_id", name)]
    df$year <- t
    rows[[as.character(t)]] <- df[, c("hex_id", "year", name)]
  }
  do.call(rbind, rows)
}


#' hist_infest : somme des sévérités TBE observées sur t-5..t-1 (par année)
#'
#' Dérivée du panel cible (sévérité 0=Aucun .. k). Pour chaque année t, on somme
#' la sévérité de l'hexagone sur les `window` années précédentes.
build.tbe.history <- function(target_panel, item, name) {
  window <- item$window %||% 5
  years <- sort(unique(target_panel$year))
  # Matrice hex × année de sévérité (0 si non observé).
  wide <- reshape.severity.wide(target_panel)

  rows <- list()
  for (t in years) {
    prev_years <- as.character((t - window):(t - 1))
    prev_years <- intersect(prev_years, colnames(wide))
    val <- if (length(prev_years) == 0) rep(0, nrow(wide))
           else rowSums(wide[, prev_years, drop = FALSE], na.rm = TRUE)
    df <- data.frame(hex_id = as.integer(rownames(wide)), year = t)
    df[[name]] <- val
    rows[[as.character(t)]] <- df
  }
  do.call(rbind, rows)
}


#' dist_foyer : distance (centroïde) au foyer actif le plus proche à t-1
#'
#' "Actif" = sévérité correspondant à `active_levels`. Pour l'année t, on repère
#' les hexagones actifs en t-1 et on calcule, pour chaque hexagone, la distance
#' de son centroïde au centroïde actif le plus proche (0 si l'hexagone est
#' lui-même actif). NA s'il n'y a aucun foyer actif en t-1.
build.tbe.distance <- function(target_panel, hex, item, name, active_severity) {
  years <- sort(unique(target_panel$year))
  wide <- reshape.severity.wide(target_panel)
  cent <- st_centroid(st_geometry(hex))
  hex_ids <- hex$hex_id

  rows <- list()
  for (t in years) {
    py <- as.character(t - 1)
    active_ids <- if (!py %in% colnames(wide)) integer(0)
                  else as.integer(rownames(wide))[wide[, py] >= active_severity]
    df <- data.frame(hex_id = hex_ids, year = t)
    if (length(active_ids) == 0) {
      df[[name]] <- NA_real_
    } else {
      act_cent <- cent[match(active_ids, hex_ids)]
      dmat <- st_distance(cent, act_cent)
      df[[name]] <- apply(dmat, 1, min)
    }
    rows[[as.character(t)]] <- df
  }
  do.call(rbind, rows)
}


#' Passer le panel cible (long) en matrice hex × année de sévérité (0 par défaut)
reshape.severity.wide <- function(target_panel) {
  hex_ids <- sort(unique(target_panel$hex_id))
  years   <- sort(unique(target_panel$year))
  wide <- matrix(0, nrow = length(hex_ids), ncol = length(years),
                 dimnames = list(hex_ids, years))
  idx_h <- match(target_panel$hex_id, hex_ids)
  idx_y <- match(target_panel$year, years)
  wide[cbind(idx_h, idx_y)] <- target_panel$severity
  wide
}


# -----------------------------------------------------------------------------
# ASSEMBLAGE DU PANEL + CACHE
# -----------------------------------------------------------------------------

#' Construire le panel hex×année complet et l'écrire en cache
#'
#' @return (invisible) la liste mise en cache : $data, $hex, $types, $meta.
build.covariate.panel <- function(config_path = "config.yaml") {
  cfg <- read_yaml(config_path)
  target_crs <- cfg$project$target_crs %||% 32198
  min_cov <- cfg$hex_grid$min_coverage %||% 0.25
  build_cfg <- cfg$covariates_build

  # --- 1. Zone + grille ---
  zone <- load.study.zone(cfg$study_zone, target_crs)
  hex <- create.hex.grid(zone, cellsize = cfg$hex_grid$cellsize, target_crs = target_crs)

  # --- 2. Cible TBE par année (sévérité) ---
  target <- st_read(cfg$target$source, quiet = TRUE)
  target[[cfg$target$value_field]] <- normalize.intensity(
    target[[cfg$target$value_field]], cfg$target$ordinal_levels)
  target <- st_transform(target, target_crs)
  target <- target[zone, ]   # découpage à la zone (cohérent avec mod_spa_temp)
  years <- sort(unique(target[[cfg$target$year_field]]))
  levels_all <- c("Aucun", cfg$target$ordinal_levels)   # 0 = Aucun

  target_rows <- list()
  for (t in years) {
    v <- aggregate.target.year(hex, target, cfg$target, t, min_cov)
    sev <- match(as.character(v), levels_all) - 1L   # 0..k
    sev[is.na(sev)] <- 0L                             # absence = 0
    target_rows[[as.character(t)]] <- data.frame(hex_id = hex$hex_id, year = t,
                                                 severity = sev)
  }
  target_panel <- do.call(rbind, target_rows)

  # --- 3. Covariables (via les constructeurs) ---
  items <- build_cfg$items
  static_tabs <- list(); temporal_tabs <- list(); types <- character(0)
  eco_data <- NULL
  needs_eco <- any(sapply(items, function(it) grepl("^ecoforest", it$builder)))
  if (needs_eco) eco_data <- load.ecoforest(cfg)

  for (nm in names(items)) {
    it <- items[[nm]]; types[nm] <- it$value_type
    message(">> covariable : ", nm, " (", it$builder, ")")
    tab <- switch(it$builder,
      ecoforest_host_prop = build.prop.hote(hex, eco_data, it, nm),
      ecoforest_age       = build.age(hex, eco_data, it, nm, min_cov),
      winter_tmin         = build.winter.tmin(hex, cfg, it, nm, years, min_cov),
      tbe_history         = build.tbe.history(target_panel, it, nm),
      tbe_distance        = build.tbe.distance(target_panel, hex, it, nm,
                              active_severity = min(match(it$active_levels, levels_all) - 1L)),
      stop("Constructeur inconnu : ", it$builder))
    if (isTRUE(it$temporal)) temporal_tabs[[nm]] <- tab else static_tabs[[nm]] <- tab
  }

  # --- 4. Fusion en panel long (hex_id × year) ---
  panel <- target_panel
  panel$target <- factor(levels_all[panel$severity + 1L], levels = levels_all, ordered = TRUE)
  for (nm in names(static_tabs))   panel <- merge(panel, static_tabs[[nm]],   by = "hex_id", all.x = TRUE)
  for (nm in names(temporal_tabs)) panel <- merge(panel, temporal_tabs[[nm]], by = c("hex_id", "year"), all.x = TRUE)
  panel <- panel[order(panel$year, panel$hex_id), ]

  # --- 5. Cache ---
  meta <- list(zone_regions = cfg$study_zone$regions, cellsize = cfg$hex_grid$cellsize,
               years = years, covariate_types = types,
               covariate_groups = lapply(items, function(it) it$groups %||% c("group1", "group2")),
               built_at = Sys.time())
  cache <- list(data = panel, hex = hex, types = types, meta = meta)

  dir.create(dirname(build_cfg$cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, build_cfg$cache)
  # Version SIG : géométrie + covariables statiques, et le panel en table.
  hex_static <- hex
  for (nm in names(static_tabs))
    hex_static <- merge(hex_static, static_tabs[[nm]], by = "hex_id", all.x = TRUE)
  st_write(hex_static, build_cfg$cache_gpkg, layer = "hex", delete_dsn = TRUE, quiet = TRUE)
  st_write(panel, build_cfg$cache_gpkg, layer = "panel", delete_layer = TRUE, quiet = TRUE)

  message("\nPanel écrit : ", nrow(panel), " lignes (", length(unique(panel$hex_id)),
          " hexagones × ", length(years), " années)")
  message("Cache : ", build_cfg$cache, " + ", build_cfg$cache_gpkg)
  print(utils::head(panel))
  invisible(cache)
}

if (sys.nframe() == 0) {
  build.covariate.panel("config.yaml")
}
