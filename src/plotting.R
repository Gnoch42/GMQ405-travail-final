library(ggplot2)
library(ggpubr)
library(tidyr)
library(tmap)

# Aide à la décision pour le choix de matrice spatiale
moran.I.comparison <- function(df) {
  ggplot(data = df, aes(x = reorder(mat, moran.I), y = moran.I)) +
    geom_segment(aes(
      x = reorder(mat, moran.I),
      xend = reorder(mat, moran.I),
      y = 0,
      yend = moran.I
    )) +
    geom_point(size = 4,
               fill = "red",
               shape = 21) +
    xlab("Matrice de pondération spatiale") +
    ylab("I de Moran") +
    coord_flip()
}

lisa.map <- function(data, field, year) {
  colors <- c("red", "blue", "lightpink", "skyblue2", "lightgray")
  
  map <- tm_shape(data) +
    tm_polygons(
      col = 'black',
      lwd = 0.5,
      fill = field,
      fill.scale = tm_scale_categorical(values = colors),
      fill.legend = tm_legend(
        title = 'Typologie',
        position = tm_pos_in('right', 'bottom'),
        frame = FALSE
      )
    ) +
    tm_scalebar(breaks = c(0, 10, 20),
                position = tm_pos_in('left', 'bottom')) +
    tm_title(year) +
    tm_layout(frame = FALSE)
  
  tmap_save(map, paste0("data/maps/lisa_", year, ".png"), width = 5, height = 4, dpi = 300)
  
  return(map)
}

lisa.map.no.legend <- function(data, field, year) {
  colors <- c("red", "blue", "lightpink", "skyblue2", "lightgray")
  
  map <- tm_shape(data) +
    tm_polygons(
      col = 'black',
      lwd = 0.5,
      fill = field,
      fill.scale = tm_scale_categorical(values = colors),
      fill.legend = tm_legend_hide()
    ) +
    tm_title(year) +
    tm_layout(frame = FALSE)
  
  tmap_save(map, paste0("data/maps/lisa_", year, ".png"), width = 5, height = 4, dpi = 300)
  
  return(map)
}


# Aide la la décision pour l'optimisation du m dans GeoCMeans
m.comparison <- function(data) {
  # Graphique avec l'inertie expliquée
  inertie <- ggplot(data) +
    geom_line(aes(x = m, y = Explained.inertia)) +
    geom_point(aes(x = m, y = Explained.inertia), color = "red") +
    labs(title = "a. Variation des données expliquées", y = "Inertie expliquée", x = "Paramètre m")
  # Graphique avec l'indice de silhouette
  silhouette <- ggplot(data) +
    geom_line(aes(x = m, y = Silhouette.index)) +
    geom_point(aes(x = m, y = Silhouette.index), color = "red") +
    labs(title = "b. Consistance des groupes", y = "Critère de silhouette floue", x = "Paramètre m")
  # Combinaison des deux graphiques dans la figure
  ggarrange(inertie, silhouette)
}

# Aide la la décision pour l'optimisation du alpha dans GeoCMeans
alpha.comparison <- function(data) {
  silouhette <- ggplot(data) +
    geom_line(aes(x = alpha, y = Silhouette.index), color = 'black') +
    geom_point(aes(x = alpha, y = Silhouette.index), color = 'red') +
    labs(x = "Alpha", y = "Indice de silhouette")
  
  incoherance <- ggplot(data) +
    geom_line(aes(x = alpha, y = spConsistency), color = 'black') +
    geom_point(aes(x = alpha, y = spConsistency), color = 'red') +
    labs(x = "Alpha", y = "Indice d'incohérence spatiale")
  
  ggarrange(silouhette, incoherance, ncol = 2, nrow = 1)
}

# Graphique de comparaison des groupes de GeoCMeans
group.mean.comparison <- function(data, group.fields) {
  
  data_long <- data %>%
    pivot_longer(
      cols = all_of(group.fields),
      names_to = "Cluster",
      values_to = "valeur"
    )
  
  ggplot(data_long, aes(x = year, y = valeur,
                        color = Cluster, group = Cluster)) +
    geom_line() +
    geom_point() +
    labs(x = "Année", y = "Niveau TBE") +
    scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "yellow"))
}
