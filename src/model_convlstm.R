# =============================================================================
# src/model_convlstm.R
# -----------------------------------------------------------------------------
# Modèle 3 du groupe 1 : ConvLSTM (réseau spatio-temporel) via `torch`.
#
# Principe : au lieu de traiter chaque hexagone isolément (Markov, RF), le
# ConvLSTM voit la CARTE ENTIÈRE d'intensité TBE comme une image, et une SÉQUENCE
# d'années comme une vidéo. Il apprend directement la dynamique de propagation
# spatiale (où et à quelle vitesse le foyer s'étend/régresse) pour prédire la
# carte de l'année suivante.
#
# Chaîne de traitement :
#   1. Rastériser la cible (et, en option, les covariables) : grille hexagonale
#      -> pile d'images régulières (une par année).
#   2. Empiler en tenseur (échantillon, temps, canaux, hauteur, largeur).
#   3. Entraîner un ConvLSTM (convolution + mémoire récurrente) à prédire la
#      carte de classes de l'année t+1 à partir des `seq_len` années précédentes.
#   4. Prédire la ou les années de test, puis reconvertir la carte prédite en
#      valeurs par hexagone (classe majoritaire des cellules de chaque hexagone).
#
# Le ConvLSTM a besoin de GRANDES cartes et de séries longues : sur une petite
# zone il surapprend. Il est donc pertinent sur les 2 MRC (voire plus), pas sur
# une seule municipalité.
#
# NB méthodologique : la validation croisée par BLOCS SPATIAUX n'est pas définie
# pour un réseau qui raisonne sur la carte entière -> convlstm.predict.full()
# renvoie NULL dans ce cas (l'évaluation l'affiche en NA) et ne participe qu'au
# TEST TEMPOREL (sa vraie vocation : prédire le futur).
#
# Dépendances : torch (backend installé via torch::install_torch()), terra, sf.
# =============================================================================

library(torch)
library(terra)
library(sf)

`%||%` <- function(a, b) if (is.null(a)) b else a


## ARCHITECTURE (torch) ------------------------------------------------------

#' Cellule ConvLSTM : un LSTM dont les portes sont des convolutions 2D
#'
#' Les quatre portes (entrée, oubli, sortie, cellule candidate) sont calculées
#' par UNE convolution sur la concaténation [entrée, état caché], ce qui préserve
#' la structure spatiale (contrairement à un LSTM classique qui aplatirait tout).
convlstm_cell <- nn_module("convlstm_cell",
  initialize = function(in_ch, hidden_ch, kernel = 3) {
    self$hidden <- hidden_ch
    self$conv <- nn_conv2d(in_ch + hidden_ch, 4 * hidden_ch, kernel,
                           padding = kernel %/% 2)
  },
  forward = function(x, h, c) {
    g  <- self$conv(torch_cat(list(x, h), dim = 2))   # concat sur les canaux
    ch <- torch_chunk(g, 4, dim = 2)
    i <- torch_sigmoid(ch[[1]]); f <- torch_sigmoid(ch[[2]])
    o <- torch_sigmoid(ch[[3]]); g_c <- torch_tanh(ch[[4]])
    c2 <- f * c + i * g_c                              # nouvelle mémoire
    list(h = o * torch_tanh(c2), c = c2)               # nouvel état caché
  }
)

#' Réseau ConvLSTM : parcourt la séquence puis classe chaque pixel de t+1
convlstm_net <- nn_module("convlstm_net",
  initialize = function(in_ch, hidden_ch, n_classes, kernel = 3) {
    self$hidden <- hidden_ch
    self$cell <- convlstm_cell(in_ch, hidden_ch, kernel)
    self$head <- nn_conv2d(hidden_ch, n_classes, 1)   # classifieur par pixel
  },
  forward = function(x) {                              # x : (b, T, C, H, W)
    d <- dim(x); b <- d[1]; Tn <- d[2]; H <- d[4]; W <- d[5]
    h <- torch_zeros(b, self$hidden, H, W)
    c <- torch_zeros(b, self$hidden, H, W)
    for (t in seq_len(Tn)) {
      out <- self$cell(x[, t, , , ], h, c)
      h <- out$h; c <- out$c
    }
    self$head(h)                                       # (b, n_classes, H, W)
  }
)


## RASTÉRISATION : paires (hexagone × année) -> images régulières -------

#' Reconstituer la sévérité par hexagone et par année depuis les paires
#'
#' Chaque paire (hex, year_t) fournit l'état en t (state_t) ET en t+1 (state_t1),
#' ce qui permet de remplir toutes les années, y compris les trous (ex. l'année
#' de validation exclue du découpage).
severity.wide <- function(all_pairs, hex) {
  years_t <- sort(unique(all_pairs$year_t))
  present <- sort(unique(c(years_t, years_t + 1)))
  full    <- seq(min(present), max(present))          # plage consécutive complète
  sev <- matrix(NA_integer_, nrow = nrow(hex), ncol = length(full),
                dimnames = list(as.character(hex$hex_id), as.character(full)))
  for (y in years_t) {
    p <- all_pairs[all_pairs$year_t == y, ]
    sev[as.character(p$hex_id), as.character(y)]     <- p$state_t
    sev[as.character(p$hex_id), as.character(y + 1)] <- p$state_t1
  }
  # Une année sans aucune donnée (ex. l'année de validation, exclue du découpage)
  # crée un trou dans la séquence : on la comble en reportant la dernière carte
  # connue (LOCF), et on remplit par 0 (absence) les hexagones non observés des
  # années présentes.
  observed <- apply(!is.na(sev), 2, any)
  last <- NULL
  for (j in seq_len(ncol(sev))) {
    if (observed[j]) { sev[is.na(sev[, j]), j] <- 0L; last <- j }
    else if (!is.null(last)) sev[, j] <- sev[, last]
  }
  first <- which(observed)[1]                          # backfill éventuel en tête
  if (first > 1) for (j in seq_len(first - 1)) sev[, j] <- sev[, first]
  storage.mode(sev) <- "integer"
  sev
}

#' Gabarit raster + correspondance cellule -> hexagone (pour aller/retour)
raster.template <- function(hex, res) {
  hexv <- terra::vect(hex)
  tmpl <- terra::rast(terra::ext(hexv), resolution = res, crs = terra::crs(hexv))
  hex_r <- terra::rasterize(hexv, tmpl, field = "hex_id")
  list(tmpl = tmpl, hid = as.vector(terra::values(hex_r)),
       H = terra::nrow(tmpl), W = terra::ncol(tmpl))
}

#' Transformer un vecteur "valeur par hexagone" en image H×W (via la corresp.)
hexvals.to.image <- function(vals_by_hex, rt, fill = 0) {
  v <- vals_by_hex[as.character(rt$hid)]
  v[is.na(v)] <- fill
  matrix(v, nrow = rt$H, ncol = rt$W, byrow = TRUE)
}


## ENTRAÎNEMENT + PRÉDICTION -----------------------------------------------

#' Prédiction ConvLSTM complète (interface utilisée par le pipeline et l'éval.)
#'
#' @param train,test Paires de transition (build.transition.pairs).
#' @param n_states   Nombre d'états.
#' @param hex        Grille hexagonale (géométrie) pour la rastérisation.
#' @param cov_names  Covariables à ajouter comme canaux (si use_covariates).
#' @param cfg        Config (section `convlstm`).
#' @return Vecteur d'états prédits aligné sur `test`, ou NULL (torch absent,
#'   échec, ou mode CV spatiale — non défini pour un réseau sur carte entière).
convlstm.predict.full <- function(train, test, n_states, hex, cov_names, cfg) {
  # CV spatiale (test au sein des années d'entraînement) : non défini -> NULL.
  if (max(test$year_t) <= max(train$year_t)) return(NULL)
  # Garde-fou torch : si le backend n'est pas prêt, on renonce proprement.
  ok <- tryCatch({ torch::torch_tensor(1); TRUE }, error = function(e) FALSE)
  if (!ok) { message("  [ConvLSTM] backend torch indisponible -> ignoré."); return(NULL) }

  cc <- cfg$convlstm %||% list()
  res     <- cc$raster_res %||% 2000
  seq_len <- cc$seq_len %||% 4
  hidden  <- cc$hidden %||% 16
  epochs  <- cc$epochs %||% 60
  lr      <- cc$lr %||% 0.01
  torch::torch_manual_seed(1)

  all_pairs <- rbind(train[, intersect(names(train), names(test))],
                     test[,  intersect(names(train), names(test))])
  # Seules les covariables NUMÉRIQUES peuvent servir de canaux au réseau.
  cov_names <- cov_names[vapply(cov_names, function(cn) is.numeric(all_pairs[[cn]]), logical(1))]
  use_cov <- isTRUE(cc$use_covariates) && length(cov_names) > 0

  # --- 1. Rastérisation ---
  rt  <- raster.template(hex, res)
  sev <- severity.wide(all_pairs, hex)
  all_years <- as.integer(colnames(sev))
  mask_img <- matrix(as.integer(!is.na(rt$hid)), nrow = rt$H, ncol = rt$W, byrow = TRUE)

  # Images de sévérité (canal principal, normalisé 0..1) par année.
  sev_img <- lapply(colnames(sev), function(y)
    hexvals.to.image(setNames(sev[, y], rownames(sev)), rt) / (n_states - 1))

  # Canaux covariables (option) : standardisés globalement, valeur/hexagone/année.
  cov_img <- NULL
  if (use_cov) {
    stats <- lapply(cov_names, function(cn) {
      v <- all_pairs[[cn]]; c(m = mean(v, na.rm = TRUE), s = stats::sd(v, na.rm = TRUE))
    })
    names(stats) <- cov_names
    cov.year <- function(y) {                      # canaux covariables de l'année y
      yy <- min(y, max(all_pairs$year_t))          # années hors plage -> plus proche
      lapply(cov_names, function(cn) {
        p <- all_pairs[all_pairs$year_t == yy, c("hex_id", cn)]
        vb <- setNames(p[[cn]], p$hex_id)
        img <- hexvals.to.image(vb, rt, fill = NA)
        img <- (img - stats[[cn]]["m"]) / (stats[[cn]]["s"] + 1e-6)  # standardisation
        img[is.na(img)] <- 0                       # NA -> moyenne (0 après standard.)
        img
      })
    }
    cov_img <- lapply(all_years, cov.year)
    names(cov_img) <- as.character(all_years)
  }

  n_channels <- 1 + if (use_cov) length(cov_names) else 0

  # Construit le tenseur d'entrée (seq_len, C, H, W) se terminant à l'année `end`.
  make.input <- function(end_year) {
    yrs <- (end_year - seq_len + 1):end_year
    arr <- array(0, dim = c(seq_len, n_channels, rt$H, rt$W))
    for (k in seq_along(yrs)) {
      yi <- match(as.character(yrs[k]), colnames(sev))
      arr[k, 1, , ] <- sev_img[[yi]]
      if (use_cov) for (j in seq_along(cov_names)) arr[k, 1 + j, , ] <- cov_img[[as.character(yrs[k])]][[j]]
    }
    arr
  }

  # --- 2. Échantillons d'entraînement : cibles = années de transition d'entraînement ---
  target_years <- sort(unique(train$year_t)) + 1          # année prédite = t+1
  target_years <- target_years[(target_years - seq_len) >= min(all_years)]
  if (length(target_years) == 0) { message("  [ConvLSTM] série trop courte -> ignoré."); return(NULL) }

  X <- array(0, dim = c(length(target_years), seq_len, n_channels, rt$H, rt$W))
  Y <- array(1L, dim = c(length(target_years), rt$H, rt$W))   # classes 1..C
  for (i in seq_along(target_years)) {
    ty <- target_years[i]
    X[i, , , , ] <- make.input(ty - 1)
    Y[i, , ] <- hexvals.to.image(setNames(sev[, as.character(ty)], rownames(sev)), rt) + 1L
  }
  x_t <- torch::torch_tensor(X, dtype = torch::torch_float())
  y_t <- torch::torch_tensor(Y, dtype = torch::torch_long())
  m_t <- torch::torch_tensor(array(rep(mask_img, length(target_years)),
                                   dim = c(length(target_years), rt$H, rt$W)),
                             dtype = torch::torch_float())

  # Pondération des classes (déséquilibre : l'absence domine largement).
  freq <- tabulate(as.vector(Y[m_t$to(dtype = torch::torch_bool()) |> as.array()]), nbins = n_states)
  w <- sum(freq) / (n_states * pmax(freq, 1)); w <- w / mean(w)
  w_t <- torch::torch_tensor(w, dtype = torch::torch_float())

  # --- 3. Entraînement ---
  model <- convlstm_net(n_channels, hidden, n_states)
  opt <- torch::optim_adam(model$parameters, lr = lr)
  for (e in seq_len(epochs)) {
    opt$zero_grad()
    logits <- model(x_t)                                   # (N, C, H, W)
    loss_map <- torch::nnf_cross_entropy(logits, y_t, weight = w_t, reduction = "none")
    loss <- (loss_map * m_t)$sum() / m_t$sum()
    loss$backward(); opt$step()
  }

  # --- 4. Prédiction des années de test + reconversion en hexagones ---
  hid <- rt$hid
  pred_state_by_hex_year <- list()
  for (yt in sort(unique(test$year_t))) {
    ty <- yt + 1
    inp <- torch::torch_tensor(array(make.input(yt), dim = c(1, seq_len, n_channels, rt$H, rt$W)),
                               dtype = torch::torch_float())
    logits <- model(inp)
    la <- as.array(logits$squeeze(1))                      # (C, H, W)
    cls_img <- apply(la, c(2, 3), which.max) - 1L          # classe 0..C-1
    cls_cells <- as.vector(t(cls_img))                     # ordre lignes (comme hid)
    # Classe majoritaire par hexagone.
    agg <- tapply(cls_cells[!is.na(hid)], hid[!is.na(hid)], function(z)
      as.integer(names(which.max(table(z)))))
    pred_state_by_hex_year[[as.character(yt)]] <- agg
  }

  # Aligner sur les lignes de `test`.
  vapply(seq_len(nrow(test)), function(i) {
    agg <- pred_state_by_hex_year[[as.character(test$year_t[i])]]
    v <- agg[as.character(test$hex_id[i])]
    if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v)
  }, integer(1))
}
