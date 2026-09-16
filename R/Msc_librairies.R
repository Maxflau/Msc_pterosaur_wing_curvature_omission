# ============================================================
# Ce script charge tous les "packages" (extensions R) nécessaires
# à la pipeline. Il faut les avoir installés au moins une fois
# avec install.packages("nom_du_package").
# ============================================================
library(dplyr)      # manipuler des tableaux de données
library(ggplot2)     # créer des graphiques
library(tidyr)       # remettre en forme des tableaux
library(magick)       # lire/traiter des images
library(parallel)     # calculs en parallèle (plus rapide)
library(readr)        # lire des fichiers texte/CSV
library(ggpubr)       # graphiques statistiques
# --- Morphométrie (analyse de forme) ---
library(Momocs)
library(sfsmisc)
library(geomorph)
library(morphospace)
library(mvMORPH)
# --- Phylogénie (arbres évolutifs) ---
library(ape)
library(phytools)
library(nlme)
# --- Analyse fonctionnelle ---
library(alphahull)
library(pracma)
library(MASS)
# --- Statistiques ---
library(dunn.test)
library(pairwiseAdonis)
library(knitr)
library(tibble)
library(randomForest)
library(pls)
library(forcats)
# --- Graphiques ---
library(gridExtra)
library(grid)
library(RColorBrewer)
library(scales)
library(viridis)
library(ggnewscale)
library(ggrepel)
library(cowplot)
library(akima)
library(mgcv)
library(magrittr)
library(patchwork)
library(viridisLite)

cat("✓ All libraries loaded\n\n")
# Charge les fonctions maison utilisées pour extraire les contours des ailes
source("R/MorphometricExtraction_Functions.r")
source("R/MorphoFiles_Function.r")