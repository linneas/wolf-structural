#################
### written by LS
### Plotting SFS for all unrelated wolves
### (using only deletions, pre genotype frequency filter and manual curation)

rm(list=ls())
d<-Sys.Date()
format(d, format="%Y-%m-%d")

################################################################################
# Setting up, loading R libraries
require(vcfR)
require(tidyverse)
require(viridis)
require(gtable)
require(grid)
library(ggh4x)

plotdir<-"plots/"


################################################################################
# Input files
meta_file="helpfiles/metadata.txt"
unrel_file="helpfiles/unrelated_ind.txt"
all_pol_file="data/polarization/joint.all.txt"
un_vcf_file="data/filtered/duphold.fold.MHSQ.DEL.50-10000.vcf"
un_vep_file="vep/DEL.duphold50-10000.comb.txt"

#Read in data
un_vcf <- read.vcfR(un_vcf_file)
meta_tib<-meta_file %>% read.table(header=TRUE) %>% as_tibble() %>% rename(Indiv=UU_ID, Cat=Category)
ind_tib<-unrel_file %>% read.table(header=FALSE) %>% as_tibble() %>% rename(Indiv=V1)
un_vep_tib <- un_vep_file %>% read.table(header=FALSE) %>% as_tibble() %>%
							rename(ID=V1, Type=V2) %>% mutate(ID=as.character(ID))
pol_tib <- all_pol_file %>% read.table(header=FALSE) %>% as_tibble() %>%
 							rename(ID=V1, AA=V2) %>% mutate(ID=as.character(ID))
# Convert to tidy format
un_tidy_vcf <- vcfR2tidy(un_vcf, format_fields=c("GT"), dot_is_NA=TRUE)


################################################################################
# USING DELETERIOUS AND INTRONS (BOTH AA=0/0 AND AA=1/1)

# All unrelated, all polarized sites
unrel_allele_cnt <- un_tidy_vcf$gt %>% filter(!is.na(gt_GT)) %>%
	inner_join(un_tidy_vcf$fix) %>% inner_join(meta_tib) %>%
	inner_join(ind_tib) %>%
	inner_join(un_vep_tib) %>% inner_join(pol_tib) %>%
	select(CHROM,POS,ID,Type,Indiv,gt_GT,AA,SVLEN,SVTYPE) %>%
	mutate(new_gt=case_when((AA=='1/1' & gt_GT=='0/0') ~'1/1',
	                      (AA=='1/1' & gt_GT=='1/1') ~'0/0',
												TRUE ~ gt_GT)) %>%
	select(ID, Indiv, new_gt, Type) %>%
	group_by(ID, new_gt, Type) %>% summarize(count=n()) %>%
	mutate(ancestral=case_when((new_gt=='0/0') ~ (count*2.0),
															(new_gt=='0/1') ~ (count*1.0), TRUE ~ 0),
				derived=case_when((new_gt=='1/1') ~ (count*2.0),
													(new_gt=='0/1') ~ (count*1.0), TRUE ~ 0)) %>%
	group_by(ID, Type) %>%
	summarize(totanc=sum(ancestral), totder=sum(derived))

# Make sfs using all polarized sites
all_unrel_sfs <- unrel_allele_cnt %>%
	filter(Type=="deleterious" | Type=="intronic") %>%
	mutate(Type=case_when(Type=="intronic" ~ "Intronic",
	 											Type=="deleterious" ~ "Deleterious",
												TRUE ~ Type)) %>%
	filter(totanc+totder==(2*length(ind_tib$Indiv)) && totder>0) %>%
	group_by(Type, totder) %>% summarize(count=n()) %>%
	mutate(frac=case_when((Type=='Intronic') ~ (count/sum(count[Type=='Intronic'])),
												(Type=='Deleterious') ~ (count/sum(count[Type=='Deleterious'])))) %>%
	ungroup()

# Binning, all polarised
all_bin_sfs<-all_unrel_sfs %>% mutate(freq=totder/80) %>%
				mutate(bin=case_when(freq<0.05 ~ 0.025,
														freq>=0.05 & freq<0.1 ~ 0.075,
														freq>=0.1 & freq<0.15 ~ 0.125,
														freq>=0.15 & freq<0.2 ~ 0.175,
														freq>=0.2 & freq<0.25 ~ 0.225,
														freq>=0.25 & freq<0.3 ~ 0.275,
														freq>=0.3 & freq<0.35 ~ 0.325,
														freq>=0.35 & freq<0.4 ~ 0.375,
														freq>=0.4 & freq<0.45 ~ 0.425,
														freq>=0.45 & freq<0.5 ~ 0.475,
														freq>=0.5 & freq<0.55 ~ 0.525,
														freq>=0.55 & freq<0.6 ~ 0.575,
														freq>=0.6 & freq<0.65 ~ 0.625,
														freq>=0.65 & freq<0.7 ~ 0.675,
														freq>=0.7 & freq<0.75 ~ 0.725,
														freq>=0.75 & freq<0.8 ~ 0.775,
														freq>=0.8 & freq<0.85 ~ 0.825,
														freq>=0.85 & freq<0.9 ~ 0.875,
														freq>=0.9 & freq<0.95 ~ 0.925,
														freq>=0.95 ~ 0.975)) %>%
			group_by(Type, bin) %>%
			summarize(totfrac=sum(frac), totcount=sum(count)))


################################################################################
# PLOT

# Unrelated, own binning
p<-ggplot(all_bin_sfs, aes(x=bin, y=totfrac, fill=Type)) +
geom_col(position="dodge", alpha = 0.8, just=0) +
scale_fill_manual(values=c("lightblue", "dodgerblue4")) +
labs(x="Derived allele frequency", y="Fraction of sites")+
theme(panel.grid.major = element_line(colour = '#f5f4e6'),
			panel.background = element_rect(fill = '#FFFDF8', colour = '#f5f4e6'),
	panel.spacing = unit(2, "lines"),
	legend.position="bottom",
	legend.title=element_blank(),
	axis.ticks = element_blank(),
	strip.background = element_rect(fill='white'),
	strip.text = element_blank())

outfile=paste(plotdir,"Figure4.",d,".png", sep="")
ggsave(outfile, plot = p, scale = 1, dpi = 600, limitsize = TRUE, width=4, height=4)





################################################################################
# STATS - COMPARE DISTRIBUTIONS WITH CHI-SQUARE GOODNESS OF FIT TEST
# Compare distribution of deleterious and synonymous

all_unrel_bin_pivot <- all_bin_sfs %>% select(Type, bin, totcount, totfrac) %>%
              pivot_wider(names_from=Type, values_from=c(totcount, totfrac)) %>%
              replace_na(list(totcount_Deleterious=0, totfrac_Deleterious=0))
chisq.test(all_unrel_bin_pivot$totcount_Deleterious, p=all_unrel_bin_pivot$totfrac_Intronic)
#data:  all_unrel_bin_pivot$totcount_Deleterious
#X-squared = 55.84, df = 19, p-value = 1.727e-05
