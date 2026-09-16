##################################################################################
# LOAD LIBRARIES
##################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(magick)
library(parallel)
library(readr)
library(ggpubr)

# Morphometric libraries
library(Momocs)
library(sfsmisc)
library(geomorph)
library(morphospace)
library(mvMORPH)

# Phylogenetic libraries
library(ape)
library(phytools)
library(nlme)

# Functional analysis libraries
library(alphahull)
library(pracma)
library(MASS)

#Statistical library
library(dunn.test)
library(pairwiseAdonis)
library(knitr)
library(readr)
library(tibble)
library(randomForest)
library(pls)
library(forcats)
library(tibble)

# Plotting
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
library(ggnewscale)
library(patchwork)
library(viridisLite)

cat("✓ All libraries loaded\n\n")
# Source functions
source("./MorphometricExtraction_Functions.r")
source("./MorphoFiles_Function.r")