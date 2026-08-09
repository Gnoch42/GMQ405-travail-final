compute.historic.stats <- function(tbe, study.zone) {
  Sup.totale.ha <- as.numeric(st_area(study.zone)) / 10000
  cat("Superficie totale du territoire d'étude :", round(Sup.totale.ha), "ha\n")
  
  ## 4.1 Découpage précis des polygones à la frontière du territoire d'étude----
  # (nécessaire car SupHaCea peut inclure une superficie hors de nos 2 MRC
  tbe.clip <- st_intersection(tbe, study.zone)
  
  ## 4.2 Recalcul de la superficie après découpage, en hectares----
  tbe.clip$Sup.ha <- as.numeric(st_area(tbe.clip)) / 10000
  
  ## 4.3 Superficie par année et par sévérité (Léger/Modéré/Grave)----
  Calculs.severite <- st_drop_geometry(tbe.clip)
  Calculs.severite <- group_by(Calculs.severite, ANNEE, Niveau)
  Calculs.severite <- summarise(Calculs.severite, Superficie_ha = sum(Sup.ha), .groups = "drop")
  
  ## 4.4 Passage en format "large" : une colonne par niveau de sévérité----
  Tableau.severite <- pivot_wider(Calculs.severite,
                                  names_from = Niveau,
                                  values_from = Superficie_ha,
                                  values_fill = 0)
  
  ## 4.5 Superficie totale touchée et % du territoire touché, par année----
  
  Calculs.touche <- group_by(tbe.clip, ANNEE)
  Calculs.touche <- summarise(Calculs.touche, geom = st_union(geom), .groups = "drop")
  
  # Calcul des champs
  Calculs.touche$Superficie_touchee_ha <- as.numeric(st_area(Calculs.touche)) / 10000
  Calculs.touche$Pct_territoire_touche <- round((Calculs.touche$Superficie_touchee_ha / Sup.totale.ha) * 100, 2)
  Calculs.touche <- st_drop_geometry(Calculs.touche)
  
  ## 4.6 Nombre de polygones de perturbation par année ----
  Calculs.nb <- st_drop_geometry(tbe.clip)
  Calculs.nb <- group_by(Calculs.nb, ANNEE)
  Calculs.nb <- summarise(Calculs.nb, Nb_polygones = n(), .groups = "drop")
  
  ## 4.7 Assemblage du tableau final (années en lignes, variables en colonnes)----
  TableauFinal <- left_join(Tableau.severite, Calculs.touche, by = "ANNEE")
  TableauFinal <- left_join(TableauFinal, Calculs.nb, by = "ANNEE")
  TableauFinal <- arrange(TableauFinal, ANNEE)
  
  ## 4.8 Export du tableau en CSV ----
  write.csv(TableauFinal, "data/output/resultats_ampleur_tbe_2014_2025.csv", row.names = FALSE)
  
  return(TableauFinal)
}
