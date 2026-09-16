
# --- Outline repair functions (shared) ----------------------------------------
close_loop <- function(c0){ n<-nrow(c0); per<-sum(sqrt(diff(c0[,1])^2+diff(c0[,2])^2))
gap<-sqrt((c0[1,1]-c0[n,1])^2+(c0[1,2]-c0[n,2])^2)
if(per>0 && gap/per>0.001) c0<-rbind(c0,c0[1,,drop=FALSE]); c0 }
resample_perimeter <- function(c0,n_out=100){ x<-c0[,1];y<-c0[,2]
seg<-sqrt(diff(x)^2+diff(y)^2); cum<-c(0,cumsum(seg)); tot<-cum[length(cum)]
if(tot==0) return(c0[rep(1,n_out),,drop=FALSE]); tg<-seq(0,tot,length.out=n_out+1)[1:n_out]
cbind(x=approx(cum,x,xout=tg,ties="ordered")$y, y=approx(cum,y,xout=tg,ties="ordered")$y) }
align_start_to_tip <- function(c0){ ti<-which.max(c0[,1]); n<-nrow(c0); c0[c(ti:n,seq_len(ti-1)),,drop=FALSE] }
normalise_centroid_size <- function(c0){ c0[,1]<-c0[,1]-mean(c0[,1]); c0[,2]<-c0[,2]-mean(c0[,2])
cs<-sqrt(sum(c0[,1]^2+c0[,2]^2)); if(cs>0) c0<-c0/cs; c0 }
prepare_outline <- function(c0){ normalise_centroid_size(align_start_to_tip(resample_perimeter(close_loop(c0),100))) }

# --- Build EFA + PCA ONCE ------------------------------------------------------
shared_names <- intersect(names(outlines_list), performance_data_clean$species)
outline_matrices <- lapply(shared_names, function(nm) prepare_outline(as.matrix(outlines_list[[nm]][,c("x","y")])))
names(outline_matrices) <- shared_names

out_obj   <- Out(outline_matrices) %>% coo_center() %>% coo_close()
efa_obj   <- efourier(out_obj, nb.h=100, norm=FALSE, start=FALSE)
pca_shape <- PCA(efa_obj)
shape_var <- (pca_shape$sdev^2/sum(pca_shape$sdev^2))*100
lab1 <- sprintf("Shape PC1 (%.1f%%)",shape_var[1]); lab2 <- sprintf("Shape PC2 (%.1f%%)",shape_var[2])
cat(sprintf("Shape PC1: %.1f%% | PC2: %.1f%%\n",shape_var[1],shape_var[2]))

shape_scores <- data.frame(species=rownames(pca_shape$x), shapePC1=pca_shape$x[,1], shapePC2=pca_shape$x[,2]) %>%
  left_join(performance_data_clean %>% dplyr::select(species,clade,aspect_ratio), by="species")

# --- Shared helpers ------------------------------------------------------------
spc1_range <- range(shape_scores$shapePC1); spc2_range <- range(shape_scores$shapePC2)
calculate_hull <- function(x,y){ if(length(x)<3) return(data.frame(x=x,y=y))
  p<-cbind(x,y); i<-chull(p); data.frame(x=p[c(i,i[1]),1],y=p[c(i,i[1]),2]) }
build_hulls <- function(){ h<-data.frame()
for(cl in unique(na.omit(shape_scores$clade))){ cd<-shape_scores[!is.na(shape_scores$clade)&shape_scores$clade==cl,]
if(nrow(cd)>=3){ z<-calculate_hull(cd$shapePC1,cd$shapePC2); z$clade<-cl; h<-rbind(h,z) } }
names(h)<-c("shapePC1","shapePC2","clade"); h }
shape_hulls <- build_hulls()

# Place real silhouettes at given coordinates (cols x_coord,y_coord per specimen).
place_silhouettes <- function(coord_df, scale){ ov<-data.frame()
for(nm in shared_names){ r<-coord_df[coord_df$species==nm,]; if(nrow(r)==0) next
c0<-outline_matrices[[nm]]; sx<-c0[,1]-mean(c0[,1]); sy<-c0[,2]-mean(c0[,2])
me<-max(abs(c(sx,sy))); if(me>0){sx<-sx/me*scale; sy<-sy/me*scale}
ov<-rbind(ov, data.frame(x=sx+r$x_coord[1], y=sy+r$y_coord[1], species=nm, clade=r$clade[1])) }; ov }

# Reconstruct a theoretical wing by inverse EFA at PC scores (pc1,pc2).
mean_coe<-pca_shape$center; rotation<-pca_shape$rotation; nb_h<-length(mean_coe)/4; n_pc<-ncol(rotation)
recon_wing <- function(pc1,pc2,nb_pts=90){ s<-numeric(n_pc); s[1]<-pc1; s[2]<-pc2
ce<-as.numeric(mean_coe+s%*%t(rotation))
m<-as.matrix(efourier_i(list(an=ce[1:nb_h],bn=ce[(nb_h+1):(2*nb_h)],cn=ce[(2*nb_h+1):(3*nb_h)],dn=ce[(3*nb_h+1):(4*nb_h)]),nb.pts=nb_pts))
data.frame(x=m[,1],y=m[,2]) }
build_grid_wings <- function(nx,ny,scale_frac){ e1<-diff(spc1_range)*0.06; e2<-diff(spc2_range)*0.06
gx<-seq(spc1_range[1]-e1,spc1_range[2]+e1,length.out=nx); gy<-seq(spc2_range[1]-e2,spc2_range[2]+e2,length.out=ny)
g<-expand.grid(shapePC1=gx,shapePC2=gy); ws<-((max(gx)-min(gx))/(nx-1))*scale_frac; out<-vector("list",nrow(g))
for(i in seq_len(nrow(g))){ w<-recon_wing(g$shapePC1[i],g$shapePC2[i]); w$x<-w$x-mean(w$x); w$y<-w$y-mean(w$y)
me<-max(abs(c(w$x,w$y))); if(me>0){w$x<-w$x/me*ws; w$y<-w$y/me*ws}
w$gx<-w$x+g$shapePC1[i]; w$gy<-w$y+g$shapePC2[i]; w$grid_id<-i; out[[i]]<-w }; do.call(rbind,out) }

# --- FIGURE 1: real silhouettes at shape coordinates --------------------------
sil_coords <- transform(shape_scores, x_coord=shapePC1, y_coord=shapePC2)
sil <- place_silhouettes(sil_coords, diff(spc1_range)/22)
p1 <- ggplot() +
  geom_polygon(data=shape_hulls, aes(shapePC1,shapePC2,group=clade), fill=NA, colour="grey70", linewidth=0.3) +
  geom_polygon(data=sil, aes(x,y,group=species,fill=clade), colour="grey30", linewidth=0.15, alpha=0.7) +
  labs(title="Shape morphospace — real silhouettes", x=lab1, y=lab2, fill="Clade") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave("output/plots_PDF/Morphospace_silhouettes.pdf", p1, width=9, height=8, dpi=300)

# --- FIGURE 2: theoretical EFA grid + specimen points -------------------------
grid_w <- build_grid_wings(8,7,0.2)
p2 <- ggplot() +
  geom_polygon(data=grid_w, aes(gx,gy,group=grid_id), fill="grey82", colour="grey50", linewidth=0.1) +
  geom_polygon(data=shape_hulls, aes(shapePC1,shapePC2,group=clade), fill=NA, colour="black", linewidth=0.15) +
  geom_point(data=shape_scores, aes(shapePC1,shapePC2,colour=clade), size=2) +
  labs(title="Theoretical shape morphospace (EFA)", x=lab1, y=lab2, colour="Clade") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave("output/plots_PDF/Morphospace_theoretical_grid.pdf", p2, width=10, height=8, dpi=300)

# --- FIGURE 3: dense Foth-style grid + open black points ----------------------
grid_dense <- build_grid_wings(18,14,0.35)
p3 <- ggplot() +
  geom_polygon(data=grid_dense, aes(gx,gy,group=grid_id), fill="grey72", colour=NA) +
  geom_point(data=shape_scores, aes(shapePC1,shapePC2), shape=21, fill="white", colour="black", size=1.5, stroke=0.5) +
  labs(title="Empirical Shape Distribution", x=lab1, y=lab2) +
  coord_equal(expand=FALSE) + theme_classic(base_size=12) +
  theme(plot.title=element_text(hjust=0.5,size=16))
ggsave("output/plots_PDF/Morphospace_Fothstyle.pdf", p3, width=9, height=8, dpi=300)

