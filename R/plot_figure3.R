#################
### Plotting PCA from plink

#Clear the workspace
rm(list=ls())

# Setup Load packages
require(tidyverse)
library(viridis)
library(ggh4x)
library(cowplot)
plotdir="plots/"
d<-Sys.Date()
format(d, format="%Y-%m-%d")


################################################################################
# Input files
meta_file="helpfiles/metadata.txt"
del_cur_vec_file="plink/curated.1cur.strict.DEL.pca.eigenvec"
dup_cur_vec_file="plink/curated.1cur.strict.DUP.pca.eigenvec"
inv_cur_vec_file="plink/curated.1cur.strict.INV.pca.eigenvec"
snp_vec_file="plink/SNPs.pca.eigenvec"

del_cur_val_file="plink/curated.1cur.strict.DEL.pca.eigenval"
dup_cur_val_file="plink/rejected.1cur.strict.DUP.pca.eigenval"
inv_cur_val_file="plink/rejected.1cur.strict.INV.pca.eigenval"
snp_val_file="plink/SNPs.pca.eigenval"

#Diagonal sums from the Relation covariance matrices, needed for percent explained
sum_file="plink/final.diagsum"


# Read in the data, start with eigenval and diagonal sums to get
# percent variance explained
del_cur_eig<-del_cur_val_file %>% read.table(header=FALSE) %>% as_tibble()
dup_cur_eig<-dup_cur_val_file %>% read.table(header=FALSE) %>% as_tibble()
inv_cur_eig<-inv_cur_val_file %>% read.table(header=FALSE) %>% as_tibble()
snp_eig<-snp_val_file %>% read.table(header=FALSE) %>% as_tibble()
sum_rel<-sum_file %>% read.table(header=FALSE, colClasses=c("character", "numeric"))

# Add PC and % explained as columns to eigenval tibble
del_cur_eig <- tibble::rowid_to_column(del_cur_eig, "PC")
del_cur_eig <- del_cur_eig %>% mutate(Expl=100*V1/sum_rel$V2[sum_rel$V1=="DEL"])
dup_cur_eig <- tibble::rowid_to_column(dup_cur_eig, "PC")
dup_cur_eig <- dup_cur_eig %>% mutate(Expl=100*V1/sum_rel$V2[sum_rel$V1=="DUP"])
inv_cur_eig <- tibble::rowid_to_column(inv_cur_eig, "PC")
inv_cur_eig <- inv_cur_eig %>% mutate(Expl=100*V1/sum_rel$V2[sum_rel$V1=="INV"])
snp_eig <- tibble::rowid_to_column(snp_eig, "PC")
snp_eig <- snp_eig %>% mutate(Expl=100*V1/sum_rel$V2[sum_rel$V1=="SNP"])
# Extract % Explained for PC1 and PC2, to be used later
p1del <- del_cur_eig %>% filter(PC==1) %>% select(Expl) %>% pull()
p2del <- del_cur_eig %>% filter(PC==2) %>% select(Expl) %>% pull()
p1dup <- dup_cur_eig %>% filter(PC==1) %>% select(Expl) %>% pull()
p2dup <- dup_cur_eig %>% filter(PC==2) %>% select(Expl) %>% pull()
p1inv <- inv_cur_eig %>% filter(PC==1) %>% select(Expl) %>% pull()
p2inv <- inv_cur_eig %>% filter(PC==2) %>% select(Expl) %>% pull()
p1snp <- snp_eig %>% filter(PC==1) %>% select(Expl) %>% pull()
p2snp <- snp_eig %>% filter(PC==2) %>% select(Expl) %>% pull()


# PC AND METADATA
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble()
del_cur_tib<-del_cur_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Deletions", Num="25,640", PCex1=p1del, PCex2=p2del)
dup_cur_tib<-dup_cur_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Duplications", PC1=PC1*-1, PC2=PC2*-1, Num="786",
              PCex1=p1dup, PCex2=p2dup)
inv_cur_tib<-inv_cur_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3) %>%
              mutate(Type="Inversions", Num="126", PCex1=p1inv, PCex2=p2inv)
snp_tib<-snp_vec_file %>% read.table(header=FALSE) %>% as_tibble() %>%
              rename(UU_ID=V1, PC1=V2, PC2=V3)  %>%
              mutate(Type="SNPs", Num="8,313,538", PCex1=p1snp, PCex2=p2snp)

# Merge and add metadata
comb_tib <- bind_rows(del_cur_tib, dup_cur_tib, inv_cur_tib, snp_tib) %>%
      select(UU_ID, PC1, PC2, Type, Num, PCex1, PCex2) %>% inner_join(meta_tib) %>%
      mutate(Population=case_when((Category=="Finland" | Category=="Russia") ~ Category,
                                  Category=="NR_immigrants" ~ "Non-reproducing immigrants",
                                  Category=="R_immigrants" ~ "Reproducing immigrants",
                                  TRUE ~ "Scandinavia")) %>%
      select(UU_ID, PC1, PC2, Type, Batch, Population, Num, PCex1, PCex2)

# Convert columns to factors
comb_tib$Population <- factor(comb_tib$Population, levels = c("Scandinavia", "Finland", "Russia", "Reproducing immigrants", "Non-reproducing immigrants"))
comb_tib$Type<-factor(comb_tib$Type, levels=c("SNPs","Deletions","Duplications","Inversions"))

# Extract Founder
found<-comb_tib %>% filter(UU_ID=="D-85-01")
################################################################################
# Plot PC1 vs PC2

# Make plot without actually plotting
p_orig<-ggplot(comb_tib, aes(x=PC1, y=PC2, shape=Population, color=Population, fill=Population)) +
  geom_point(alpha=0.5) +
  facet_wrap2(vars(Type), nrow = 2, ncol = 2, trim_blank = FALSE) +
  #labeller = labeller(Type = list(PC1 = "PC1 Label", PC2 = "PC2 Label"))
  scale_color_viridis(discrete=TRUE) +
	scale_fill_viridis(discrete=TRUE) +
	scale_shape_manual(values=c(21, 22, 23, 24, 25))+
  geom_point(data=found, aes(x=PC1, y=PC2), shape=21, color="red", fill="red", show.legend=FALSE, alpha=0.7)+
  geom_text(aes(x=-0.18, y=0.32, label=Num), size=3, hjust=0, color="black")+
  geom_text(aes(x=0, y=-0.30, label=paste("(",round(PCex1,1),"%)", sep="")), size=3, hjust=0.5, color="darkgrey")+
  geom_text(aes(x=-0.18, y=0, label=paste("(",round(PCex2,1),"%)", sep="")), size=3, hjust=0,vjust=0.5, color="darkgrey")+
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
        panel.spacing = unit(1, "lines"),
				axis.ticks.x=element_blank(),
				axis.ticks.y=element_blank(),
        legend.position = "bottom",
				legend.key = element_rect(fill = "#FFFDF8"),
        legend.title = element_blank()) +
  guides(color=guide_legend(nrow=3))


#####
# To get separate axis labels (with % explained), Make four different plots and
# combine. Take legend from standard facet plot above

myplot <- function(x) {
  comb_tib %>%
    filter(Type == x) %>%
    ggplot(aes(x=PC1, y=PC2, shape=Population, color=Population, fill=Population)) +
      geom_point(alpha=0.5) +
      facet_wrap(vars(Type), scales = "free") +
      scale_color_viridis(discrete=TRUE) +
    	scale_fill_viridis(discrete=TRUE) +
    	scale_shape_manual(values=c(21, 22, 23, 24, 25))+
      geom_point(data=found[which(found$Type==x),], aes(x=PC1, y=PC2), shape=21, color="red", fill="red", show.legend=FALSE, alpha=0.7)+
      theme(panel.grid.major = element_line(colour = '#f5f4e6'),
            panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
            strip.text.x = element_text(size = 12),
    				axis.ticks.x=element_blank(),
    				axis.ticks.y=element_blank(),
            legend.position = "none")+
      labs(x = paste("PC1 (",round(exp_df$PC1[exp_df$Type==x],1),"%)", sep=""), y = paste("PC2 (",round(exp_df$PC2[exp_df$Type==x],1),"%)", sep=""))
  }

# Get legend from a joint plot above
legend <- get_legend(
p_orig +
  guides(color = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = "#FFFDF8"),
        legend.title = element_blank()))

#Add the figures generated with the function into a single plot
p_new<-lapply(unique(comb_tib$Type), myplot) %>%
  plot_grid(plotlist = ., nrow = 2, ncol=2)

#Organize the plots together with the legend (legend is 20% of plot height)
p_fin<-plot_grid(p_new, legend, ncol = 1, rel_heights = c(1, .2))+theme(plot.background = element_rect(fill = "white"))

outfile=paste(plotdir,"Figure3.",d,".png", sep="")
ggsave(outfile, plot = p_fin, scale = 1, dpi = 600, limitsize = TRUE,
      width=6, height=6.5)
