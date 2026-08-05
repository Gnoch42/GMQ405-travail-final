library(ggplot2)

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