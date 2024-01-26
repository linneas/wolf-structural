#################
### written by LS
### Plotting PCA for Kept, Removed (both quality filtered and manually rejected)
### and Uncurated (did not pass genotype frequency filter)

#Clear the workspace
rm(list=ls())

# Setup Load packages
require(tidyverse)
library(viridis)
plotdir="plots/"
d<-Sys.Date()
format(d, format="%Y-%m-%d")


################################################################################
# Using 1 curator

# Input files and read in the data
meta_file="helpfiles/metadata.txt"
prefix=paste(cur,filt,sep=".")
del_rem_vec_file="plink/rejected_removed.DEL.pca.eigenvec"
del_unc_vec_file="plink/uncurated.DEL.pca.eigenvec"
del_cur_vec_file="plink/curated.1cur.strict.DEL.pca.eigenvec"
dup_rem_vec_file="plink/rejected_removed.DUP.pca.eigenvec"
dup_unc_vec_file="plink/uncurated.DUP.pca.eigenvec"
dup_cur_vec_file="plink/curated.1cur.strict.DEL.pca.eigenvec"
inv_rem_vec_file="plink/rejected_removed.INV.pca.eigenvec"
inv_unc_vec_file="plink/uncurated.INV.pca.eigenvec"
inv_cur_vec_file="plink/curated.1cur.strict.INV.pca.eigenvec"

meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% mutate(Batch=as.character(Batch))
del_rem_tib<-del_rem_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Deletions", Set="Rejected")
del_unc_tib<-del_unc_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Deletions", Set="Uncurated")
del_cur_tib<-del_cur_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Deletions", Set="Kept")
dup_rem_tib<-dup_rem_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Set="Rejected")
dup_unc_tib<-dup_unc_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Set="Uncurated")
dup_cur_tib<-dup_cur_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Set="Kept")
inv_rem_tib<-inv_rem_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Inversions", Set="Rejected")
inv_unc_tib<-inv_unc_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3)  %>%
              mutate(Type="Inversions", Set="Uncurated")
inv_cur_tib<-inv_cur_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3)  %>%
              mutate(Type="Inversions", Set="Kept")

# Merge and add metadata
comb_tib <- bind_rows(del_cur_tib, del_rem_tib, del_unc_tib, dup_cur_tib, dup_rem_tib,
          dup_unc_tib, inv_cur_tib, inv_rem_tib, inv_unc_tib) %>%
          select(UU_ID, PC1, PC2, Type, Set) %>%
          inner_join(meta_tib) %>% mutate(Population=if_else(Category=="Finland" | Category=="Russia" | Category=="NR_immigrants" |Category=="R_immigrants", Category, "Scandinavia")) %>%
          select(UU_ID, PC1, PC2, Type, Set, Population, Batch)

# Convert columns to factors
comb_tib$Population <- factor(comb_tib$Population, levels = c("Scandinavia", "Finland", "Russia", "R_immigrants", "NR_immigrants"))
comb_tib$Set<-factor(comb_tib$Set, levels=c("Kept","Rejected","Uncurated"))
comb_tib$Type<-factor(comb_tib$Type, levels=c("Deletions","Duplications","Inversions"))
comb_tib$Batch<-factor(comb_tib$Batch, levels=c("2016","2017","2019"))


################################################################################
# Plot PC1 vs PC2

p<-ggplot(comb_tib, aes(x=PC1, y=PC2, shape=Population, color=Batch, fill=Batch)) +
  geom_point(alpha=0.5) +
   ggh4x::facet_grid2(Set~Type, scales="free",independent="x") +
  scale_color_viridis(discrete=TRUE) +
	scale_fill_viridis(discrete=TRUE) +
	scale_shape_manual(values=c(21, 22, 23, 24, 25))+
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
				axis.ticks.x=element_blank(),
				axis.ticks.y=element_blank(),
				legend.key = element_rect(fill = "#FFFDF8"))

outfile=paste(plotdir,"FigureS2.",d,".jpg", sep="")
ggsave(outfile, plot = p,scale = 1,dpi = 300,limitsize = TRUE,
  width=8,height=7)
