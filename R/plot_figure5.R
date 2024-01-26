#################
### written by LS
### Plotting load per generation


################################################################################
# Setting up, loading R libraries
require(vcfR)
require(tidyverse)
library(ggpubr)
require(cowplot)
d<-Sys.Date()
format(d, format="%Y-%m-%d")
plotdir<-"plots/"

################################################################################
# Input files
meta_file="helpfiles/metadata.txt"
del_vcf_file="data/curated/1cur.strict.DEL.vcf"
del_vep_file="vep/DEL.comb.detail.txt"
dup_vcf_file="data/curated/1cur.strict.DUP.vcf"
dup_vep_file="vep/DUP.comb.detail.txt"
inv_vcf_file="data/curated/1cur.strict.INV.vcf"
inv_vep_file="vep/INV.comb.detail.txt"
all_pol_file="data/polarization/joint.all.txt"


#Read in data
del_vcf <- read.vcfR(del_vcf_file)
dup_vcf <- read.vcfR(dup_vcf_file)
inv_vcf <- read.vcfR(inv_vcf_file)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
del_vep_tib <- del_vep_file %>% read.table(header=FALSE) %>% as_tibble() %>%
							rename(ID=V1, Type=V2) %>% mutate(ID=as.character(ID))
dup_vep_tib <- dup_vep_file %>% read.table(header=FALSE) %>% as_tibble() %>%
							rename(ID=V1, Type=V2) %>% mutate(ID=as.character(ID))
inv_vep_tib <- inv_vep_file %>% read.table(header=FALSE) %>% as_tibble() %>%
							rename(ID=V1, Type=V2) %>% mutate(ID=as.character(ID))

all_pol_tib <- all_pol_file %>% read.table(header=FALSE) %>% as_tibble() %>%
 							rename(ID=V1, AA=V2) %>% mutate(ID=as.character(ID))

# Convert to tidy format
del_tidy_vcf <- vcfR2tidy(del_vcf, format_fields=c("GT"), dot_is_NA=TRUE)
dup_tidy_vcf <- vcfR2tidy(dup_vcf, format_fields=c("GT"), dot_is_NA=TRUE)
inv_tidy_vcf <- vcfR2tidy(inv_vcf, format_fields=c("GT"), dot_is_NA=TRUE)

# Merge the vcfR data into one big tibble
del_all<-del_tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	       inner_join(del_tidy_vcf$fix) %>% inner_join(del_vep_tib)
dup_all<-dup_tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	       inner_join(dup_tidy_vcf$fix) %>% inner_join(dup_vep_tib)
inv_all<-inv_tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	       inner_join(inv_tidy_vcf$fix) %>% inner_join(inv_vep_tib)
sv_all<-bind_rows(del_all,dup_all,inv_all) %>% inner_join(meta_tib) %>%
			inner_join(all_pol_tib)


# Extracting relevant data and combine tibbles
gt_tib <- sv_all %>% filter(Type!="noncoding") %>%
	mutate(Type=case_when((Type=="full_overlap" | Type=="partial_CDS") ~ "Deleterious",
												Type=="intronic" ~ "Intronic",
												TRUE ~ Type)) %>%
	select(CHROM,POS,ID,Type,Indiv,gt_GT,AA,SVLEN,SVTYPE) %>%
	mutate(new_gt=case_when((AA=='1/1' & gt_GT=='0/0') ~'1/1',
                          (AA=='1/1' & gt_GT=='1/1') ~'0/0',
													TRUE ~ gt_GT)) %>%
	group_by(Indiv,Type, AA, new_gt) %>% summarize(count=n()) %>% ungroup() %>%
	mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
	pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
	replace_na(list(gt00=0, gt01=0, gt11=0)) %>%
	mutate(sum=gt00+gt01+gt11) %>%
	mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, homancfrq=gt00/sum, deralfreq=(gt11*2+gt01)/(sum*2))

extr_gt <-gt_tib %>% inner_join(meta_tib) %>%
	select(Indiv,Type,Cat,hetfrq,homderfrq,deralfreq,homancfrq, AA,GenClass) %>%
	mutate(Grp=case_when(GenClass=='Finland' | GenClass=='Russia' | GenClass=='NR_immigrants' ~ 'Reference',
				TRUE ~ 'Scandinavia')) %>% filter(AA=="00" & GenClass!="Founders")

# Make factor
extr_gt$GenClass <- factor(extr_gt$GenClass, levels = c("Finland", "Russia", "NR_immigrants", "F1", "F2", "F3", "F4", "F5", "F6", "R_immigrants", "L1", "L2", "L3"))
extr_gt$Type <- factor(extr_gt$Type, levels= c("Deleterious", "Intronic"))

################################ PLOTTING ######################################

# a) Heterozygous sites
pA<-ggplot(extr_gt,aes(x=GenClass, y=hetfrq))+
geom_point(data=extr_gt, colour="#440154", size=1, alpha=0.3) +
stat_summary(fun=mean, geom="point", colour="#440154", shape="-", size=10, alpha=0.9) +
facet_grid(Type~factor(Grp, level=c("Reference","Scandinavia")), scale="free_x", space="free_x") +
  labs(x=NULL,y="Proportion of heterozygous derived genotypes")+
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust=1),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12))

# b) Homozygous derived sites
pB<-ggplot(extr_gt,aes(x=GenClass, y=homderfrq))+
geom_point(data=extr_gt, colour="#440154", size=1, alpha=0.3) +
stat_summary(fun=mean, geom="point", colour="#440154", shape="-", size=10, alpha=0.9) +
facet_grid(Type~factor(Grp, level=c("Reference","Scandinavia")), scales="free_x", space="free_x") +
  labs(x=NULL,y="Proportion of homozygous derived genotypes")+
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust=1),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12))


# Combine a) and b)
pcomb <- plot_grid(pA,pB, nrow = 2, ncol=1, labels = c("(a)", "(b)"))

output <- paste(plotdir,"Figure5.",d,".png", sep="")
ggsave(output, plot = pcomb, scale = 1, unit = "in", dpi = 600, limitsize = TRUE,
	width = 7.5, height = 10)


################################################################################
# Numbers for text:
extr_gt %>% select(homderfrq, Type, GenClass) %>% group_by(Type,GenClass) %>%
summarize(mean=mean(homderfrq), sd=sd(homderfrq)) %>% print(n=50)
