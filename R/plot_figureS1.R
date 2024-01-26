#################
### written by LS
### Plotting PCA for Duplications and Inversions, comparing 1 and 2 curators

#Clear the workspace
rm(list=ls())

# Setup, oad packages
require(tidyverse)
library(viridis)
library(ggh4x)
plotdir="plots/"
d<-Sys.Date()
format(d, format="%Y-%m-%d")


################################################################################
# Input files
meta_file="helpfiles/metadata.txt"
dup_cur1_vec_file="plink/curated.1cur.strict.DUP.pca.eigenvec"
inv_cur1_vec_file="plink/curated.1cur.strict.INV.pca.eigenvec"
dup_cur2_vec_file="plink/curated.2cur.strict.DUP.pca.eigenvec"
inv_cur2_vec_file="plink/curated.2cur.strict.INV.pca.eigenvec"
dup_rej1_vec_file="plink/rejected.1cur.strict.DUP.pca.eigenvec"
inv_rej1_vec_file="plink/rejected.1cur.strict.INV.pca.eigenvec"
dup_rej2_vec_file="plink/rejected.2cur.strict.DUP.pca.eigenvec"
inv_rej2_vec_file="plink/rejected.2cur.strict.INV.pca.eigenvec"

# Read in the data
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble()
dup_cur1_tib<-dup_cur1_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Cur="One curator", Set="Kept", PC1=PC1*-1, Num="786")
inv_cur1_tib<-inv_cur1_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Inversions",Cur="One curator", Set="Kept",  Num="126")
dup_cur2_tib<-dup_cur2_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Cur="Two curators", Set="Kept", PC1=PC1*-1, Num="413")
inv_cur2_tib<-inv_cur2_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Inversions",Cur="Two curators", Set="Kept",  Num="71")
dup_rej1_tib<-dup_rej1_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Cur="One curator", Set="Rejected", PC1=PC1*-1, Num="1496")
inv_rej1_tib<-inv_rej1_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Inversions",Cur="One curator", Set="Rejected",  Num="424")
dup_rej2_tib<-dup_rej2_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", Cur="Two curators", Set="Rejected", PC1=PC1*-1, Num="1869")
inv_rej2_tib<-inv_rej2_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Inversions",Cur="Two curators", Set="Rejected", PC2=PC2*-1, Num="479")


# Merge and add metadata
comb_tib <- bind_rows(dup_cur1_tib, inv_cur1_tib, dup_cur2_tib, inv_cur2_tib,
  dup_rej1_tib, inv_rej1_tib, dup_rej2_tib, inv_rej2_tib) %>%
      select(UU_ID, PC1, PC2, Type, Num, Cur, Set) %>% inner_join(meta_tib) %>%
      mutate(Population=case_when((Category=="Finland" | Category=="Russia") ~ Category,
                                  Category=="NR_immigrants" ~ "NR immigrants",
                                  Category=="R_immigrants" ~ "R immigrants",
                                  TRUE ~ "Scandinavia")) %>%
      select(UU_ID, PC1, PC2, Type, Batch, Population, Num, Cur, Set)

# Convert columns to factors
comb_tib$Population <- factor(comb_tib$Population, levels = c("Scandinavia", "Finland", "Russia", "R immigrants", "NR immigrants"))
comb_tib$Batch<-factor(comb_tib$Batch, levels=c("2016","2017","2019"))
comb_tib$Type<-factor(comb_tib$Type, levels=c("Duplications","Inversions"))
comb_tib$Cur<-factor(comb_tib$Cur, levels=c("One curator","Two curators"))
comb_tib$Set<-factor(comb_tib$Set, levels=c("Kept","Rejected"))

################################################################################
# Plot PC1 vs PC2, ONLY KEPT VARIANTS
extr_tib<-comb_tib %>% filter(Set!="Rejected")

p<-ggplot(extr_tib, aes(x=PC1, y=PC2, shape=Population, color=Population, fill=Population)) +
  geom_point(alpha=0.5) +
  facet_grid(Cur~Type) +
  scale_color_viridis(discrete=TRUE) +
	scale_fill_viridis(discrete=TRUE) +
	scale_shape_manual(values=c(21, 22, 23, 24, 25))+
  geom_text(aes(x=-0.12, y=0.30, label=Num), size=3, hjust=0, color="black")+
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
				axis.ticks.x=element_blank(),
				axis.ticks.y=element_blank(),
        legend.position = "right",
				legend.key = element_rect(fill = "#FFFDF8"))

outfile=paste(plotdir,"FigureS1.",d,".jpg", sep="")
ggsave(outfile, plot = p, scale = 1, dpi = 300, limitsize = TRUE,
      width=6, height=4.5)
