#########################################################
# Test de création de modèles prédictifs spatiotemporels#
#########################################################     
rm(list = ls())
graphics.off()

# Importations
library("sf")
library("yaml")

cfg <- read_yaml("config.yaml")


# Test de modélisation de l'évolution spatiale
if(!file.exists(cfg$extraction$TBE$outpath) | cfg$extraction$TBE$force_extraction){
  source("src/data_extraction.R")
  extract.tbe.data(cfg)
}

# Test de Markov


# Random Forest/XGBoost

# ConvLSTM
