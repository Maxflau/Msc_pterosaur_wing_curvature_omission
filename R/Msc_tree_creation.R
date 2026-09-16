# ============================================================
# Construit un arbre phylogénétique complet à partir de l'arbre
# de référence "Henry", en y ajoutant les espèces manquantes.
# Chaque espèce manquante est ajoutée selon 1 de ces 4 méthodes :
#    (A) ORTHOGRAPHIE      -> renomme une pointe déjà existante
#    (B) MÊME GENRE        -> placée à côté d'une espèce du même genre
#    (C) ARBRES DE RÉFÉRENCE -> placée d'après Pegas ou Andres (précis)
#    (D) CLADE DE SECOURS  -> utilisé si l'espèce est absente de tout
#  Ne pas modifier les noms d'espèces ci-dessous : ce sont des
#  correspondances taxonomiques précises établies manuellement.
# ============================================================
library(phytools)
tr <- read.nexus("data/phylogenetics/Henry")
if (inherits(tr, "multiPhylo")) {
  tr <- tr[[1]]                       # several MPTs in file -> take the first
  attr(tr, "TipLabel") <- NULL        # drop shared-label compression metadata
  class(tr) <- "phylo"
}
tr <- reorder(tr)                     # clean internal indexing
if (any(duplicated(tr$tip.label)))    # safety: collapse accidental dup labels
  stop("Henry tree itself has duplicated tip labels: ",
       paste(unique(tr$tip.label[duplicated(tr$tip.label)]), collapse=", "))
if (!is.binary(tr)) tr <- multi2di(tr)

# ---- (A) ORTHOGRAPHIC FIXES -------------------------------
rename_map <- c(
  "Gallodactylus_canjuersensis"="Cycnorhamphus_canjuersensis",     # generic synonym, already in tree
  "Hatzegopteryx_type_only"="Hatzegopteryx_thambema",              # placeholder tip = H. thambema
  "Shenzhoupterus_sanyainus"="Meilifeilong_sanyainus",
  "Altmuehlopterus_rhamphastinus"="Germanodactylus_rhamphastinus",
  "Jidapterus_edentus_"="Jidapterus_edentus",
  "Haopterus_gracilis_"="Haopterus_gracilis",
  "Hamipterus_tianshanensis_"="Hamipterus_tianshanensis",
  "Plataleorhynchus_streptophorodo"="Plataleorhynchus_streptophorodon",
  "Gladocephaloideus_jingangshanen"="Gladocephaloideus_jingangshanensis")    # data taxon = tree tip (genus reassigned)
for (old in names(rename_map)) {
  i <- match(old, tr$tip.label); if (!is.na(i)) tr$tip.label[i] <- rename_map[[old]]
}

# ---- (B) CONGENERIC SISTERS -------------------------------
congeneric <- character(0)   # (Gladocephaloideus truncation corrected upstream)

# ---- (C) REFERENCE-INFORMED (Pegas / Andres) --------------
#Place sister taxon
ref_sister <- c(
  "Simurghia_lamegoi"        = "Simurghia_robusta",                 
  "Tupuxuara_leonardii"      = "Tupuxuara_longicristatus",          
  "Thalassodromeus_oberlii"  = "Thalassodromeus_sethi",             
  "Haliskia_peterseni"       = "Ferrodraco_lentoni",                
  "Anhanguera_santanae"      = "Anhanguera_blittersdorffi",        
  "Anhanguera_robustus"      = "Anhanguera_piscator",              
  "Liaoxipterus_brachyognathus" = "Istiodactylus_latidens",         
  "Eosipterus_yangi"         = "Beipiaopterus_chenianus",           
  "Kepodactylus_insperatus"  = "Ardeadactylus_longicollum",         
  "Tsogtopteryx_mongoliensis"   = "Cryodrakon_boreas",              # user: before Cryodrakon
  "Gobiazhdarcho_tsogtbaatari"  = "Nipponopterus_mifunensis",       # user: sister of Nipponopterus
  "Garudapterus_buffetauti"     = "Lusognathus_almadrava",          # user: sister of Lusognathus
  "Longchengopterus_zhaoi"      = "Nurhachius_ignaciobritoi",       # user: sister of Nurhachius ignaciobritoi
  "Eotephradactylus_mcintirae"  = "Seazzadactylus_venieri",        
  "Pachagnathus_benitoi"        = "Raeticodactylus_filisurensis",   
  "Makrodactylus_oligodontus"        = "Skiphosoura_bavarica",
  "Yixianopterus_jingangshanensis" = "Mimodactylus_libanensis")    


#Implement a phylogenetic branch
ref_mrca <- list(
  Gladocephaloideus_jingangshanensis = c("Feilongus_youngi","Moganopterus_zhuiana"),      # pegas (tree spelling)
  Kryptodrakon_progenitor = c("Changchengopterus_pani","Propterodactylus_frankerlae"),    # user: between these two
  Balaenognathus_maeuseri = c("Gnathosaurus_subulatus","Plataleorhynchus_streptophorodon"),# pegas (tree spelling, truncated)
  Orientognathus_chaoyangensis = c("Angustinaripterus_longicephalus","Dorygnathus_banthensis"),
  Puntanipterus_globosus = c("Noripterus_parvus","Noripterus_complicidens"),
  Gladocephaloideus_jingangshanensis = c("Feilongus_youngi","Moganopterus_zhuiana")) # andres
# Gladocephaloideus has a real congeneric tip, prefer that -> drop ref_mrca entry
ref_mrca[["Gladocephaloideus_jingangshanensis"]] <- NULL

# ---- (D) CLADE-MRCA FALLBACK (absent from all refs) -------
clade_anchor <- list(
  Rhamphinion_jenkinsi        = c("Dimorphodon_macronyx","Caelestiventus_hanseni"),
  Cascocauda_rong             = c("Anurognathus_ammoni","Vesperopterylus_lamadongensis"),
  Ceoptera_evansae            = c("Darwinopterus_modularis","Kunpengopterus_sinensis"),
  Archaeoistiodactylus_linglongtaensis = c("Darwinopterus_modularis","Kunpengopterus_sinensis"),
  Laueropterus_vitriolus  = c("Skiphosoura_bavarica","Makrodactylus_oligodontus"),
  Bergamodactylus_wildi  = c("MCSNB_8950","Carniadactylus_rosenfeldi"),
  Jianchangopterus_zhaoianus  = c("Darwinopterus_modularis","Kunpengopterus_sinensis"))

# ---- GRAFTING HELPERS -------------------------------------
graft_sister <- function(tree, new_tip, ref_tip) {
  i <- match(ref_tip, tree$tip.label)
  if (is.na(i)) { warning("missing ref tip: ", ref_tip); return(tree) }
  el <- tree$edge.length[which(tree$edge[,2]==i)]
  bind.tip(tree, new_tip, where=i, position=el/2)
}
graft_mrca <- function(tree, new_tip, anchors) {
  anchors <- anchors[anchors %in% tree$tip.label]
  if (length(anchors) < 2) { warning("anchors missing for ", new_tip); return(tree) }
  bind.tip(tree, new_tip, where=getMRCA(tree, anchors))
}

for (nt in names(congeneric)) tr <- graft_sister(tr, nt, congeneric[[nt]])
for (nt in names(ref_sister)) tr <- graft_sister(tr, nt, ref_sister[[nt]])
for (nt in names(ref_mrca))   tr <- graft_mrca(tr, nt, ref_mrca[[nt]])
for (nt in names(clade_anchor)) tr <- graft_mrca(tr, nt, clade_anchor[[nt]])

# ---- CHECK & EXPORT ---------------------------------------
added <- c(names(congeneric), names(ref_sister), names(ref_mrca), names(clade_anchor))
cat("Tips after grafting :", Ntip(tr), "\n")
cat("Successfully added  :", sum(added %in% tr$tip.label), "/", length(added), "\n")
fail <- added[!(added %in% tr$tip.label)]
if (length(fail)) cat("FAILED:\n", paste(" -", fail, collapse="\n"), "\n")
class(tr) <- "phylo"                  # ensure single-tree class before writing
attr(tr, "TipLabel") <- NULL
write.tree(tr, file="output/phylogenetics/Henry_updated.tre") 
write.nexus(tr, file="output/phylogenetics/Henry_updated.nex")
cat("Saved -> Henry_updated.nex (and .tre)\n")