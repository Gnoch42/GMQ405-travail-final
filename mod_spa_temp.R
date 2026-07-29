#########################################################
# Test de création de modèles prédictifs spatiotemporels ----
#########################################################     
rm(list = ls())
graphics.off()

# Importations
library("sf")
library("yaml")

cfg <- read_yaml("config.yaml")

source("src/data_extraction.R")
source("src/hex_grid.R")

extract.tbe.data(cfg)

data_2mrc <- st_read("data/output/tbe_data_2mrc.gpkg")

# Test de modélisation de l'évolution spatiale

# Test de Markov


# Random Forest/XGBoost

# ConvLSTM
