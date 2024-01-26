#################
### written by LS
### Plotting load for descendants to original founders compared to immigrant
### descendants to immigrants from the same time period


################################################################################
# Setting up, loading R libraries and set working directory
require(vcfR)
require(tidyverse)
library(ggpubr)
plotdir<-"plots/"
rm(list=ls())
d<-Sys.Date()
format(d, format="%Y-%m-%d")


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
	select(CHROM,POS,ID,Type,Indiv,gt_GT,AA,SVLEN,SVTYPE, Cat) %>%
	mutate(new_gt=case_when((AA=='1/1' & gt_GT=='0/0') ~'1/1',
                          (AA=='1/1' & gt_GT=='1/1') ~'0/0',
													TRUE ~ gt_GT)) %>%
	group_by(Indiv,Type, AA, new_gt) %>% summarize(count=n()) %>% ungroup() %>%
	mutate_if(is.character, str_replace_all, pattern = "/", replacement = "") %>%
	pivot_wider(names_from=new_gt, names_prefix="gt", values_from=count) %>%
	replace_na(list(gt00=0, gt01=0, gt11=0)) %>%
	mutate(sum=gt00+gt01+gt11) %>%
	mutate(hetfrq=gt01/sum, homderfrq=gt11/sum, homancfrq=gt00/sum, deralfreq=(gt11*2+gt01)/(sum*2))

	################################################################################
# Extract relevant data: Only Before and After Immigration!

# Use 2007-2014S and 2007-2014I, but remove offspring to Tiveden (=Finnish)
extr_gt <-gt_tib %>% inner_join(meta_tib) %>%
	filter(Indiv!="104-G100-14" & Indiv!="89-G67-15") %>%
	select(Indiv,Type,Cat,hetfrq,homderfrq,deralfreq,homancfrq,AA,GenClass) %>%
	mutate(Grp=case_when(Cat=='2007-2014S'  ~ 'Before',
											Cat == '2007-2014I' ~ 'After',
											TRUE ~ NA)) %>%
	filter(!is.na(Grp)) %>% filter(AA=="00") %>% #filter(GenClass!="F3") %>%
	pivot_longer(cols=c(hetfrq,homderfrq), names_to="Freqtype", values_to="Freq") %>%
	mutate(Freqtype=ifelse(Freqtype=="hetfrq", "Masked load", "Realized load"))


# Make factor
extr_gt$Grp <- factor(extr_gt$Grp, levels = c("Before", "After"))
extr_gt$Type <- factor(extr_gt$Type, levels= c( "Deleterious", "Intronic"))
extr_gt$Freqtype <- factor(extr_gt$Freqtype, levels=c("Masked load", "Realized load"))

########################### STATISTICAL TEST ###################################
# Statistical Test
# start with heterozygous, putatively neutral
# Shapiro-Wilk normality test for Without
with(extr_gt, shapiro.test(Freq[Grp == "Without" & Type=="Intronic" & Freqtype=="Masked load"]))# p = 0.1
# Shapiro-Wilk normality test for With
with(extr_gt, shapiro.test(Freq[Grp == "With" & Type=="Intronic" & Freqtype=="Masked load"]))# p = 0.1
# Both are greater than 0.05, not significantly different from normaldistr
res.ftest <- var.test(Freq[Type=="Intronic" & Freqtype=="Masked load"] ~ Grp[Type=="Intronic" & Freqtype=="Masked load"], data = extr_gt)
# p>0.05, no significant difference in variance
# => do normal t-test
set1<-extr_gt %>% filter(Grp == "Without" & Type=="Intronic" & Freqtype=="Masked load") %>% select(Freq)
set2<-extr_gt %>% filter(Grp == "With" & Type=="Intronic" & Freqtype=="Masked load") %>% select(Freq)
res <- t.test(set1, set2, var.equal = TRUE)
# Significantly different!
#or
res2 <-wilcox.test(set1$Freq, set2$Freq, alternative = "two.sided")
# Also significant!
# BUT - We don't have to run this manually for every facet, we can just add
# stat_compare_means() to the ggplot command below!!


################################ PLOTTING ######################################
# Compare means and display the significance level
# (** means <0.01, *** means <0.0001)

pboth<- ggplot(extr_gt, aes(x=Grp, y=Freq, color=Grp,fill=Grp))+
	geom_boxplot(alpha=0.5, show.legend=FALSE) +
  geom_point(data=extr_gt, size=1, alpha=0.5,show.legend=FALSE) +
	stat_compare_means(aes(label = ..p.signif..),label.x=1.5, label.y=0.35, vjust=1, show.legend=FALSE) +
	scale_fill_manual(values=c("#440154", "#73D055")) +
	scale_color_manual(values=c("#440154", "#73D055")) +
  facet_grid(Type~Freqtype) +
  labs(x=NULL,y="Proportion of genotypes")+
  theme(panel.grid.major = element_line(colour = '#f5f4e6'),
        panel.background = element_rect(fill = '#FFFDF8', colour="black"), #colour = '#f5f4e6'),
        axis.text.x = element_text(vjust=0.5, margin=margin(t=0)),
				strip.text.x = element_text(size=12),
				strip.text.y = element_text(size=12))

output <- paste(plotdir,"Figure6.",d,".png", sep="")
ggsave(output,plot = pboth,scale = 1,unit = "in",dpi = 600,limitsize = TRUE,
		width = 7.5, height = 6)
