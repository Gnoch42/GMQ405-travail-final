#########################################################
# Test de création de modèles prédictifs spatiotemporels#
#########################################################     
rm(list = ls())
graphics.off()

# Importations
library("sf")
library("yaml")

cfg <- read_yaml("config.yaml")

source("src/data_extraction.R")
source("hex_grid.R")


# Test de modélisation de l'évolution spatiale

# Test de Markov


# Random Forest/XGBoost

# ConvLSTM
