library(dplyr)

aggregate_tbe_levels_to_hex <- function(study.zone, data, cfg) {
  source("src/hex_grid.R")
  grid <- create.hex.grid(study.zone, cfg$hex_grid$cellsize)
  
  for (year in seq(cfg$lisa_analysis$start, cfg$lisa_analysis$end, 1)) {
    field.name <- paste0(cfg$target$value_field, "_",year)
    grid <- aggregate.to.hex(
      hex = grid,
      source = subset(data, ANNEE == year),
      value_field = cfg$target$value_field,
      value_type = "categorical",
      out_name = field.name
    )
  }
  return(grid)
}

convert_levels_col <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "0",
    x %in% c("Leger", "Léger") ~ "1",
    x %in% c("Modere", "Modéré") ~ "2",
    x == "Grave" ~ "3",
    .default = x
  ) |>
    as.integer()
}

convert_level <- function(x) {
  x <- x %>%
    dplyr::mutate(
      dplyr::across(starts_with("Niveau_"), convert_levels_col)
    )
  return(x)
}
