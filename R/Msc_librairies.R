# ============================================================
# This script loads all the R "packages" (extensions) needed by
# the pipeline. They must be installed at least once with
# install.packages("package_name").
# ============================================================
library(dplyr)      # handle data tables
library(ggplot2)     # create plots
library(tidyr)       # reshape data tables
library(magick)       # read/process images
library(parallel)     # run calculations in parallel (faster)
library(readr)        # read text/CSV files
library(ggpubr)       # statistical plots
# --- Morphometrics (shape analysis) ---
library(Momocs)
library(sfsmisc)
library(geomorph)
library(morphospace)
library(mvMORPH)
# --- Phylogenetics (evolutionary trees) ---
library(ape)
library(phytools)
library(nlme)
# --- Functional analysis ---
library(alphahull)
library(pracma)
library(MASS)
# --- Statistics ---
library(dunn.test)
library(pairwiseAdonis)
library(knitr)
library(tibble)
library(randomForest)
library(pls)
library(forcats)
# --- Plotting ---
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
# Load the custom functions used to extract wing outlines
source("R/MorphometricExtraction_Functions.r")
source("R/MorphoFiles_Function.r")