library(ggplot2)
library(ggpubr)
library(tidyr)
library(tmap)
library(gifski)

map.tbe <- function(tbe, study.zone) {
  tmap_mode("plot")
  
  CarteTBE2025 <-
    tm_shape(study.zone) +
    tm_polygons(fill = "grey95",
                col = "black",
                lwd = 1) +
    tm_shape(tbe) +
    tm_polygons(
      col = NA,
      fill = "Niveau",
      fill.scale = tm_scale_categorical(
        values = c(
          "Léger" = "#ffffb2",
          "Modéré" = "#fd8d3c",
          "Grave" = "#bd0026"
        )
      ),
      fill.legend = tm_legend(
        title = "Intensité de la défoliation",
        frame = FALSE,
        position = c(0.7, 0.2)
      )
    ) +
    tm_layout(frame = FALSE) +
    tm_scalebar(position = tm_pos_in("right", "bottom"),
                breaks = c(0, 20),
                text.size = 0.7) +
    tm_credits(
      "Source : MRNF, Données Québec\nAuteur : Jean Yves, Antoine Fortier, Rémy Billette",
      position = tm_pos_out("center", "bottom")
    )
  
  return(CarteTBE2025)
}

tbe.animation <- function(tbe, study.zone) {
  tbe.clip <- st_intersection(tbe, study.zone)
  
  CarteAnimee <-
    tm_shape(study.zone) +
    tm_polygons(fill = "grey95", col = "black") +
    tm_shape(tbe.clip) +
    tm_polygons(
      col = NA,
      fill = "Niveau",
      fill.scale = tm_scale_categorical(
        values = c("Léger" = "#ffffb2", "Modéré" = "#fd8d3c", "Grave" = "#bd0026")
      ),
      fill.legend = tm_legend(title = "Intensité", frame = FALSE)
    ) +
    tm_animate(frames = "ANNEE", fps = 1L, play = "loop") +
    tm_title("Évolution de l'épidémie de TBE") +
    tm_layout(frame = FALSE)
  
  tmap_animation(CarteAnimee, filename = "data/maps/animation_tbe.gif")
  
  return(CarteAnimee)
}

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
  colors <- c("red", "lightpink", "blue", "skyblue2", "lightgray")
  
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
    tm_scalebar(breaks = c(0, 20),
                text.size = 1.2,
                position = tm_pos_in("RIGHT", "BOTTOM")) +
    tm_title(year) +
    tm_layout(frame = FALSE)
  
  tmap_save(map, paste0("data/maps/lisa_", year, "legend.png"), width = 5, height = 4, dpi = 300)
  
  return(map)
}

lisa.map.no.legend <- function(data, field, year) {
  colors <- c("red", "lightpink", "blue", "skyblue2", "lightgray")
  
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
