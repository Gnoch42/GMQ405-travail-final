library(spdep)
library(sf)
library(dplyr)

lisa_analysis <- function(data, W, cfg) {
  # Centrage reduction
  for (year in seq(cfg$lisa_analysis$start, cfg$lisa_analysis$end, 1)) {
    in.field.name <- paste0(cfg$target$value_field, "_", year)
    
    # Calcul des variables centrées réduite et spatialement décalées
    z.field.name <- paste0("z_", year)
    wz.field.name <- paste0("wz_", year)
    data[[z.field.name]] <- scale(st_drop_geometry(data[[in.field.name]]))[, 1]
    data[[wz.field.name]] <- lag.listw(W, data[[z.field.name]])
    
    # Calcul des I de moran
    typo.field.name <- paste0("typo_", year)
    localMoranI <- localmoran(data[[in.field.name]], W)
    alpha <- cfg$lisa_analysis$alpha
    
    typologie <- attributes(localMoranI)$quadr$mean
    typologie <- case_when(localMoranI[, 5] < alpha ~ typologie, TRUE ~ 'Non sign.')
    
    data[[typo.field.name]] <- typologie
    data[[typo.field.name]] <- factor(
      data[[typo.field.name]],
      levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Non sign."),
      labels = c("HH", "HL", "LL", "LH", "Non sign.")
    )
  }
  
  return(data)
}

#' Test les différentes matrices de pondération spatiale et renvoie un Data.Frame des résultats
#' @param data : Data.Frame sf avec un seul champ de valeur à tester
test_W_mat <- function(data) {
  stopifnot(length(st_drop_geometry(data)) == 1)
  
  # Contiguïté
  W.Rook <- poly2nb(data, queen = FALSE) |>
    nb2listw(zero.policy = TRUE)
  W.Queen <- poly2nb(data) |>
    nb2listw(zero.policy = TRUE)
  
  # Distance
  distances <- data |>
    st_centroid() |>
    st_coordinates() |>
    dist(method = 'euclidean') |>
    as.matrix()
  diag(distances) <- 0
  
  W.InvDist <- ifelse(distances == 0, 0, 1 / distances) |>
    mat2listw(style = 'W')
  
  W.InvDist2 <- ifelse(distances == 0, 0, 1 / distances^2) |>
    mat2listw(style = 'W')
  
  # Test des matrices
  mat.name.list <- c('Rook', 'Queen', 'InvDist', 'InvDist2')
  mat.list <- list(W.Rook, W.Queen, W.InvDist, W.InvDist2)
  moran.Is <- c()
  p.values <- c()
  
  for (i in 1:length(mat.list)) {
    test <- moran.mc(as.vector(st_drop_geometry(data)[, 1]), mat.list[[i]], nsim = 999)
    moran.Is[i] <- test$statistic
    p.values[i] <- test$p.value
  }
  
  ## Construction du DF
  moran.df <- data.frame(mat = mat.name.list,
                         moran.I = moran.Is,
                         p.value = p.values)
  return(moran.df)
}

compute.group.means <- function(data, grouping.fields, mean.fields, cfg) {
  output <- data.frame(matrix(ncol = length(grouping.fields), nrow = length(mean.fields)), row.names = mean.fields)
  names(output) <- grouping.fields
  output["year"] <- as.numeric(gsub("\\D", "", mean.fields))
  
  for (field in mean.fields) {
    for (group in grouping.fields) {
      wheights <- ifelse(data[[group]] > cfg$geocmeans$group_mean_thr, data[[group]], 0)
      mean <- weighted.mean(data[[field]], wheights)
      output[field, group] <- mean
    }
  }
  
  return(output)
}
