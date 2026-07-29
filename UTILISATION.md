# Guide d'utilisation — Pipeline de modélisation spatio-temporelle TBE

Ce document explique **comment configurer, lancer et interpréter** le pipeline
de modélisation de la tordeuse des bourgeons de l'épinette (TBE). Il est écrit
pour un public **non spécialiste de l'analyse spatiale** : les termes techniques
sont définis au fil du texte et regroupés dans un [glossaire](#glossaire) à la fin.

> Le `README.md` décrit les **données sources** à télécharger. Le présent guide
> décrit **le fonctionnement du code**. Les deux sont complémentaires.

---

## 1. Vue d'ensemble

Le projet répond à une question : **où sera la TBE l'année prochaine, et
pourquoi ?** Pour y répondre, le code fait deux choses très différentes,
organisées en **deux groupes de modèles comparés séparément** (on ne compare
jamais un modèle du groupe 1 à un modèle du groupe 2 : ils ne font pas le même
travail et n'utilisent pas les mêmes métriques).

| | **Groupe 1 — Prédiction** | **Groupe 2 — Compréhension** |
|---|---|---|
| Question | Où sera la TBE en **t+1** ? | Quelles covariables sont liées à l'intensité **à un instant t** ? |
| Dimension | Spatio-**temporelle** | Coupe transversale (**une année**) |
| Modèles | Markov spatial, Random Forest/XGBoost, ConvLSTM | MCO, SLX, SAR, SEM, Durbin spatial, Durbin erreur, GAM |
| Métriques | kappa pondéré, exactitude, RMSE sur rangs | AIC/BIC, pseudo-R², test de Moran |

Une **grille hexagonale** commune sert d'unité spatiale : toutes les données
(intensité TBE, covariables) y sont ramenées, quelle que soit leur forme
d'origine (polygones ou rasters).

### Architecture des fichiers

```
config.yaml               Tous les paramètres (le SEUL fichier à éditer au quotidien)
mod_spa_temp.R            Programme principal : orchestre les 2 groupes de modèles
src/
├── hex_grid.R            Grille hexagonale + agrégation spatiale générique
├── features.R            Assemblage cible + covariables sur la grille
├── models_group1.R       Modèles de prédiction t+1
├── models_group2.R       Modèles d'analyse à un instant t
├── covariate_evaluation.R  OUTIL AUTONOME de sélection des covariables
├── simulate_data.R       Générateur de données factices (démo)
└── data_extraction.R     Extraction des données réelles (existant)
```

---

## 2. Prérequis (packages R)

Indispensables (le pipeline tourne avec ceux-ci) :

```r
install.packages(c("sf", "terra", "yaml", "spdep", "spatialreg", "mgcv",
                   "ggplot2", "corrplot"))
```

Optionnels (activent des modèles supplémentaires du groupe 1) :

```r
install.packages(c("ranger", "xgboost"))   # Random Forest / XGBoost
install.packages("keras3")                  # ConvLSTM (deep learning)
```

Si un package optionnel est absent, le modèle correspondant est **ignoré
proprement** (le pipeline ne plante pas ; le tableau de résultats indique
`package_manquant`).

---

## 3. Démarrage rapide (données simulées)

Le pipeline est **prêt à l'emploi sans aucune donnée réelle** : par défaut, il
génère un jeu factice réaliste (un foyer d'infestation qui se déplace et grossit
d'année en année) pour que tout soit testable immédiatement.

```bash
Rscript mod_spa_temp.R
```

Vous obtenez deux tableaux de comparaison, un par groupe, écrits dans
`data/output/models/` :

- `comparaison_groupe1.csv`
- `comparaison_groupe2.csv`

Et pour tester l'outil de sélection des covariables :

```bash
Rscript src/covariate_evaluation.R
```

qui écrit un classement et des graphiques dans `data/output/covariate_eval/`.

---

## 4. Configuration (`config.yaml`)

**Tout se règle dans `config.yaml`. En usage normal, on ne touche pas au code.**
Copiez `config_example.yaml` en `config.yaml` et adaptez-le.

### 4.1 Passer des données simulées aux données réelles

```yaml
simulate:
  use_simulated: TRUE     # <- mettre FALSE pour utiliser vos vraies données
```

Quand `use_simulated: FALSE`, le pipeline lit la cible réelle décrite ici :

```yaml
target:
  source: "data/output/tbe_data_2mrc.gpkg"   # fichier de la cible TBE
  year_field: "ANNEE"                          # colonne de l'année
  value_field: "Niveau"                        # colonne de l'intensité
  value_type: "ordinal"
  ordinal_levels: ["Léger", "Modéré", "Grave"] # ordre CROISSANT d'intensité
```

> Les libellés d'intensité mal accentués dans la source (« Leger », « Modere »)
> sont automatiquement harmonisés vers les `ordinal_levels`.

### 4.2 Ajouter des covariables réelles (sans toucher au code)

C'est **la seule étape** pour intégrer une nouvelle covariable : ajouter une
entrée dans la section `covariates:`. Chaque covariable peut être un **raster**
ou un **polygone**, continu ou catégoriel.

```yaml
covariates:
  temperature:                       # nom libre (sera le nom de la variable)
    type: "raster"                   # raster | vector
    value_type: "continuous"         # continuous | categorical | ordinal
    path: "data/covariates/temp_hiver.tif"
    band: 1                          # (raster) numéro de bande
    agg_method: "median"

  essence:
    type: "vector"                   # polygones
    value_type: "categorical"
    path: "data/covariates/essences.gpkg"
    field: "ESSENCE"                 # (vecteur) colonne contenant la valeur
    agg_method: "mode"
```

Règles d'agrégation appliquées automatiquement :

- **Variable continue** (température, altitude…) → **médiane** pondérée.
- **Variable qualitative** (essence, classe…) → **mode** pondéré (catégorie
  dominante).

Dans les deux cas, la pondération se fait par la **surface d'intersection**
(polygones) ou la **fraction de cellule couverte** (rasters) — jamais par un
simple décompte : ce qui occupe réellement le plus de place dans l'hexagone pèse
le plus.

### 4.3 Résolution de la grille et bordures

```yaml
project:
  target_crs: 32198        # CRS projeté (mètres). 32198 = NAD83/Québec Lambert.

hex_grid:
  cellsize: 10000          # taille de cellule EN MÈTRES (10000 = hexagones ~10 km)
  min_coverage: 0.25       # couverture minimale (0-1) pour qu'un hexagone reçoive une valeur
```

- `cellsize` est en **mètres** car le pipeline reprojette d'abord vers un CRS
  projeté (`target_crs`). C'est indispensable : les données TBE sont en degrés
  (NAD83 géographique), et « 10000 degrés » n'aurait aucun sens.
- `min_coverage` gère les **hexagones de bordure** : un hexagone couvert à moins
  de 25 % par les données reçoit `NA` (valeur jugée peu fiable) au lieu d'une
  valeur calculée sur un échantillon minuscule.

### 4.4 Réglages des modèles

```yaml
models:
  group1:
    run: ["markov", "rf", "convlstm"]      # modèles à exécuter
    temporal_split:                        # découpage TEMPOREL (jamais aléatoire)
      train_years: [2014, 2020]            # apprentissage
      valid_years: [2021, 2022]            # réglage
      test_years:  [2023, 2024]            # évaluation finale
  group2:
    year_t: 2022                           # année analysée en coupe transversale
    run: ["ols", "slx", "sar", "sem", "durbin", "durbin_error", "gam"]
```

> **Pourquoi un découpage temporel et pas aléatoire ?** On veut prédire le
> futur. Si on mélangeait au hasard les années entre apprentissage et test, le
> modèle « verrait » des bouts du futur pendant l'entraînement et sa performance
> serait artificiellement gonflée. On entraîne donc sur les premières années et
> on teste sur les suivantes.

---

## 5. Choisir ses covariables : l'outil d'évaluation

**Avant** d'ajouter des covariables au modèle, utilisez l'outil autonome pour
juger lesquelles valent la peine. Il est **indépendant** du pipeline principal.

```bash
Rscript src/covariate_evaluation.R
```

Pour l'utiliser sur vos propres covariables, appelez `evaluate.covariates()`
depuis un script R en lui passant votre grille, votre cible et votre liste de
covariables (voir la fonction `demo.covariate.evaluation()` comme modèle).

**Ce qu'il produit** (`data/output/covariate_eval/`) :

- `classement_covariables.csv` : chaque covariable avec sa **force
  d'association** (0 = aucun lien, 1 = lien fort) avec l'intensité TBE, le test
  utilisé, la p-value et une interprétation (négligeable / faible / modérée /
  forte).
- `classement_covariables.png` : le même classement en barres.
- `colinearite_covariables.png` : carte de **colinéarité** (redondance) entre
  covariables continues.

**Comment lire les résultats :**

- Une covariable en **haut du classement** est associée à l'intensité TBE :
  bonne candidate.
- Une covariable **négligeable** (comme la variable « bruit » de la démo, qui
  est du hasard pur) n'apporte rien : à écarter.
- Deux covariables **fortement colinéaires** (corrélation > 0,7, ou VIF > 5)
  portent la même information : garder l'une, écarter l'autre pour ne pas
  brouiller l'interprétation des modèles.

Les tests choisis automatiquement selon le type de variable :

- **continue** vs intensité → corrélation de **Spearman** + test de
  **Kruskal-Wallis** ;
- **catégorielle** vs intensité → **V de Cramér** (dérivé du khi-deux).

(Définitions dans le [glossaire](#glossaire).)

---

## 6. Interpréter les sorties des modèles

### 6.1 Groupe 1 — prédiction t+1 (`comparaison_groupe1.csv`)

| Colonne | Signification | Lecture |
|---|---|---|
| `accuracy` | Part d'hexagones dont la classe prédite est exactement bonne | Plus haut = mieux |
| `kappa_pondere` | Accord au-delà du hasard, **pénalisant plus les grosses erreurs** (prédire « Grave » là où c'est « Aucun ») que les petites | **Métrique de référence.** 1 = parfait, 0 = hasard |
| `rmse_rang` | Erreur quadratique moyenne sur les rangs d'intensité | Plus bas = mieux |
| `n_test` | Nombre de cas dans le jeu de test | — |
| `statut` | `ok`, `package_manquant` ou `non_implemente` | — |

Le **kappa pondéré** est privilégié car l'intensité est **ordinale** : se
tromper d'une classe voisine est moins grave que de se tromper de deux classes,
et cette métrique en tient compte.

### 6.2 Groupe 2 — analyse à un instant t (`comparaison_groupe2.csv`)

| Colonne | Signification | Lecture |
|---|---|---|
| `AIC` / `BIC` | Compromis qualité d'ajustement / complexité | **Plus bas = meilleur** modèle |
| `pseudo_R2` | Part de variance expliquée | Plus haut = mieux (indicatif) |
| `moran_p_residus` | p-value du test de Moran sur les résidus du MCO | Voir ci-dessous |
| `rho` | Force de la dépendance spatiale sur la cible (SAR/Durbin) | Proche de 0 = peu de dépendance |
| `lambda` | Force de la dépendance spatiale via les erreurs (SEM/Durbin erreur) | idem |

**Le test de Moran, point clé :** il mesure l'**autocorrélation spatiale**
(tendance des lieux proches à se ressembler) dans les résidus du modèle simple
(MCO). Si `moran_p_residus` est **très petit** (≪ 0,05), il reste de la structure
spatiale non captée → cela **justifie d'utiliser les modèles spatiaux** (SAR,
SEM, Durbin) plutôt que le MCO seul. C'est exactement ce qu'on observe sur les
données simulées, fortement structurées dans l'espace.

Les modèles vont **du plus simple au plus riche** :
MCO (ignore l'espace) → SLX (ajoute les covariables des voisins) → SAR (la valeur
dépend des voisins) → SEM (l'espace agit via les erreurs) → Durbin (combinaisons
les plus générales). Comparez leurs AIC/BIC pour choisir le meilleur compromis.

---

## 7. Passer en production (résumé)

1. Télécharger les données réelles (voir `README.md`).
2. Copier `config_example.yaml` → `config.yaml`.
3. Renseigner `target:` (chemin, colonnes) et mettre `simulate.use_simulated: FALSE`.
4. Évaluer les covariables candidates avec `src/covariate_evaluation.R`.
5. Déclarer les covariables retenues dans `covariates:` du YAML.
6. Régler `hex_grid.cellsize`, `min_coverage`, et les découpages temporels.
7. Lancer `Rscript mod_spa_temp.R` et lire les deux CSV de comparaison.

---

## Glossaire

- **Agrégation spatiale** : résumer en une seule valeur, par hexagone, des
  données de forme quelconque qui le recouvrent.
- **AIC / BIC** : critères comparant des modèles ; un score **plus bas** indique
  un meilleur compromis entre qualité d'ajustement et simplicité.
- **Autocorrélation spatiale** : tendance de lieux proches à avoir des valeurs
  semblables. Mesurée par le **test de Moran**.
- **CRS (système de coordonnées)** : façon de situer les données sur la Terre.
  Un CRS **projeté** (ex. Lambert québécois, EPSG:32198) s'exprime en mètres ;
  un CRS **géographique** (ex. NAD83, EPSG:4269) en degrés.
- **Grille hexagonale** : pavage régulier en hexagones servant d'unité spatiale
  commune. Avantage sur le carré : tous les voisins sont à égale distance.
- **Kappa pondéré** : mesure d'accord entre prédiction et réalité, corrigée du
  hasard, qui pénalise davantage les grosses erreurs de classe.
- **Matrice de contiguïté / voisinage** : tableau indiquant quels hexagones sont
  voisins (ici : ceux qui se touchent). Base des modèles spatiaux.
- **Médiane** : valeur du milieu d'une série ordonnée ; robuste aux valeurs
  extrêmes (préférée à la moyenne pour l'agrégation continue).
- **Mode** : valeur (ou catégorie) la plus fréquente ; utilisé pour les
  variables qualitatives.
- **Ordinal** : variable catégorielle **ordonnée** (Aucun < Léger < Modéré <
  Grave).
- **Raster** : donnée en grille de cellules (comme une image), ex. un modèle de
  température.
- **RMSE** : racine de l'erreur quadratique moyenne ; mesure d'écart, plus basse
  = meilleure.
- **Spearman (corrélation de rangs)** : mesure si deux variables évoluent dans
  le même sens, en comparant leurs rangs. Entre −1 et +1.
- **Kruskal-Wallis** : test vérifiant si une variable continue diffère selon les
  classes d'une variable catégorielle.
- **V de Cramér** : mesure d'association (0 à 1) entre deux variables
  catégorielles.
- **VIF (facteur d'inflation de la variance)** : détecte la redondance entre
  covariables ; VIF > 5 signale une colinéarité problématique.
- **Vecteur (données vectorielles)** : données en points/lignes/**polygones**
  (ex. contours d'infestation).
```
