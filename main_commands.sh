# COMMANDS FOR STRUCTURAL VARIATION IN FENNOSCANDIAN VOLVES
# written by Linnéa Smeds and Lars S.A. Huson


#### Software versions used:
#smoove/0.2.8
#vcftools/0.1.16
#bcftools/1.17
#snakemake/7.18.2
#python/3.9.5
#pysam/0.17.0-python3.9.5
#samplot/1.3.0
#plotcritic/1.0.1
#plink2/2.00-alpha-3.7-20221024
#GATK/3.8-0
#htslib/1.12
#snakemake/5.30.1 (for first polarization step)
#vep/99
#R/4.2.1

############################## SMOVE AND DUPHOLD ###############################
#Code for running smoove (using a sigularity container):
bash/smoove.sh

# Filter to keep only autosomes
vcftools --gzvcf data/smoove/annotated/cohort.smoove.square.anno.vcf.gz --out data/smoove/annotated/cohort.smoove.square.anno.filtered --recode --recode-INFO-all --chr 1 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 2 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 3 --chr 30 --chr 31 --chr 32 --chr 33 --chr 34 --chr 35 --chr 36 --chr 37 --chr 38 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9
# After filtering, kept 154090 out of a 167664 Sites

# Remove wiltype-only sites (require at least 1 individual with an alternative allele)
bcftools view -i 'GT="0/1" || GT="1/1"' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.smoove.square.anno.filtered.recode.vcf > data/smoove/annotated/cohort.anno.filtered.GT.vcf

# Separate SV types
bcftools view -i 'SVTYPE="DEL"' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf > data/smoove/annotated/cohort.anno.nofilter.del.vcf
bcftools view -i 'SVTYPE="DUP"' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf > data/smoove/annotated/cohort.anno.nofilter.dup.vcf
bcftools view -i 'SVTYPE="INV"' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf > data/smoove/annotated/cohort.anno.nofilter.inv.vcf

### Duphold filter
mkdir data/filtered/
# deletions (DHFFC<0.7):
bcftools view -i 'SVTYPE="DEL" & DHFFC<0.7' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf >data/filtered/duphold.fold.DEL.vcf
# 61476 DEL
# duplications (DHFFC>1.3):
bcftools view -i 'SVTYPE="DUP" & DHFFC>1.3' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf >data/filtered/duphold.fold.DUP.vcf
# 6645 DUP
# inversions (no filter):
cp data/smoove/annotated/cohort.anno.filtered.GT.vcf data/filtered/duphold.fold.INV.vcf
# 9077 INV

# ADD MSHQ Filter >=3 (As in Wold et al. 2023)
for i in "DEL" "DUP" "INV"
do
  echo $i
  bcftools view -i 'MSHQ>=3 & (GT="0/1" || GT="1/1")' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/filtered/duphold.fold.$i.vcf >data/filtered/duphold.fold.MHSQ.$i.vcf
  # Count sites
  grep -v "#" data/filtered/duphold.fold.MHSQ.$i.vcf |wc
done
# 57424 DEL, 4372 DUP, 6208 INV left

# Check if there are any sites with only fixed individuals
for i in "DEL" "DUP" "INV"
do
  bcftools view -i 'MSHQ==-1' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/filtered/duphold.fold.$i.vcf >data/filtered/duphold.fixed.$i.vcf
  grep -v "#" data/filtered/duphold.fixed.$i.vcf |wc
done
# Yes: 1924 DEL, 65 DUP and 57 INV
# A quick glance gives that some are fixed 1/1 in all wolves (= different from
# dog reference), but most seems to be have mosly 0/0 individuals and just a few
# or one homozygous alternative.

# Extract removed sites from fold change / MSHQ filtering -> "removed"
for t in "DEL" "DUP" "INV"
do
  # Make list of all variants
  bcftools view -i 'SVTYPE="'$t'"' smoove/annotated/cohort.anno.filtered.GT.vcf |grep -v "#" |cut -f3 |cut -f1 -d";" |sort > data/original/$t.id.txt
  # Take inverted intersect to find variants NOT kept after filtering
  join -v1 data/original/$t.id.txt <(cut -f2 data/filtered/duphold.fold.MHSQ.$t.list |cut -f1 -d";" |sort) >data/filtered/duphold.removed.$t.id.txt
  # Extract from vcf file
  vcftools --vcf data/smoove/annotated/cohort.anno.filtered.GT.vcf --snps data/filtered/duphold.removed.$t.id.txt --recode --recode-INFO-all --out data/filtered/duphold.removed.$t
  mv data/filtered/duphold.removed.$t.recode.vcf data/filtered/duphold.removed.$t.vcf
done

# Get lengths files of raw and after duphold filtration (for supplementary table)
mkdir data/lengths
for type in "DEL" "DUP" "INV"
do
 grep -v "#" data/filtered/duphold.fold.MHSQ.$type.vcf |cut -f8 |cut -f2 -d";" |sed 's/SVLEN=//' |sed 's/-//' >data/lengths/duphold.fold.MHSQ.$type.lengths
 bcftools view -i 'SVTYPE="'$type'"' data/smoove/annotated/cohort.anno.filtered.GT.vcf |grep -v "#" |cut -f8 |cut -f2 -d";" |sed 's/SVLEN=//' |sed 's/-//' >data/lengths/raw_calls.$type.lengths
done

# Extract a set with only duphold filtered deletions to use for SFS (can't
# ute the genotypefilter/curation and filter out rare variants)
# Save only variants >50bp and <10,000bp (the length classes with fewest false positives according to the manual curation)
bcftools view -i 'SVLEN>=-10000 & SVLEN<-50' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/filtered/duphold.fold.MHSQ.DEL.vcf >data/filtered/duphold.fold.MHSQ.DEL.50-10000.vcf



########################### SAMPLOT AND PLOTCRITIC #############################
# (Deletions are used as example)
# The variants are divided by length to simplify the manual curation

# Preparation, create multiple scripts for each batch to speed up runtime
bash/prepare_samplot_scripts_del.sh

# Samplot was installed in a conda environment
conda create --name samplotcritic
source conda_init.sh
conda activate samplotcritic
conda install -c bioconda samplot

# Generate samplots for each length batch
for set in "0-50" "50-100" "100-150" "150-175" "175-200" "200-250" "250-300" "300-400" "400-500" "500-750" "750-1000" "from1000"
do
  for file in $(ls batchscripts/parts/samplot_dup_cmds.$set.*.sh)
  do
    ls $file
    sbatch $file
  done
done

# In the same conda environment, install plotcritic:
conda install -c bioconda plotcritic
### generate plotcritic websites (one for each length batch)
mkdir plotcritic
for set in "0-50" "50-100" "100-150" "150-175" "175-200" "200-250" "250-300" "300-400" "400-500" "500-750" "750-1000" "from1000"
do
  plotcritic -p plotcritic/plcr_DEL_$set -i samplot/imgs_DEL_$set/ -s -q "Is this a deletion?" -A "z":"Yes" "m":"Maybe" "c":"No"
done

# The plotcritic htmls are opened in a web browser and all samplot figures are
# assessed and judged manually by one or two curators. Reports are generated as
# tsv files.

# The following combinations were tested:
# 1) DEL, DUP, INV curated by 1 curator, keeping all Yes (called "strict" curation)
# 2) DEL, DUP, INV curated by 1 curator, keeping all Yes/Maybe (called "relaxed" curation)
# 3) DUP, INV curated by 2 curators, keeping all Yes/Yes
# 4) DUP, INV curated by 2 curators, keeping all Yes/Maybe

# For the final results in the paper, only option 1) is presented.

# Get curated sets
mkdir data/curated
f="strict"
for t in "DEL" "DUP" "INV"
do
  echo $t
  # Get a list of all pos and ids in the input vcf
  grep -v "#" data/filtered/duphold.fold.MHSQ.$t.vcf |cut -f1,2,3,8 |awk '{split($4,t,";"); print $1":"$2"-"t[3]"\t"$3}' |sed 's/END=//' >data/filtered/duphold.fold.MHSQ.$t.list
  # Get the desired sites and ids (first strict, then relaxed)
  grep $'100.0\t0.0\t0.0' plotcritic/plcr_${t}_summary_report_Linnea.tsv |awk '{s=$3+1; print $2"\t"s"\t"$4}' >plotcritic/plcr_${t}_1cur_$f.pos.txt
  awk '{print $1":"$2"-"$3}' plotcritic/plcr_${t}_1cur_$f.pos.txt |sort |join - <(sort data/filtered/duphold.fold.MHSQ.$t.list) |cut -f2 -d" " |sort -g >plotcritic/plcr_${t}_1cur_$f.id.txt
  # Extract from vcf
  vcftools --vcf data/filtered/duphold.fold.MHSQ.$t.vcf --snps plotcritic/plcr_${t}_1cur_${f}.id.txt --recode --recode-INFO-all --out data/curated/1cur.$f.$t
    mv data/curated/1cur.$f.$t.recode.vcf data/curated/1cur.$f.$t.vcf
done

# Also extracted rejected SVs for comparison
mkdir data/rejected
for t in "DEL" "DUP" "INV"
do
  echo $t
  # Get the rejected sites and ids 
  grep -e $'0.0\t100.0\t0.0' -e $'0.0\t0.0\t100.0' plotcritic/plcr_${t}_summary_report_Linnea.tsv |awk '{s=$3+1; print $2"\t"s"\t"$4}' >plotcritic/rejected_${t}_1cur_$f.pos.txt
  awk '{print $1":"$2"-"$3}' plotcritic/rejected_${t}_1cur_$f.pos.txt |sort |join - <(sort data/filtered/duphold.fold.MHSQ.$t.list) |cut -f2 -d" " |sort -g >plotcritic/rejected_${t}_1cur_$f.id.txt
  vcftools --vcf data/filtered/duphold.fold.MHSQ.$t.vcf --snps plotcritic/rejected_${t}_1cur_${f}.id.txt --recode --recode-INFO-all --out data/rejected/1cur.$f.$t
  mv data/rejected/1cur.$f.$t.recode.vcf data/rejected/1cur.$f.$t.vcf
done

# And uncurated (sites with too few individuals per genotype)
for t in "DEL" "DUP" "INV"
do
  f="strict"
  c="1cur"
  cat plotcritic/rejected_${t}_${c}_${f}.id.txt plotcritic/plcr_${t}_${c}_${f}.id.txt |sort |join -v2 -1 1 -2 2 - <(sort -k2,2 data/filtered/duphold.fold.MHSQ.$t.list) >data/uncurated/$c.$f.$t.id.txt
  vcftools --vcf data/filtered/duphold.fold.MHSQ.$t.vcf --snps data/uncurated/$c.$f.$t.id.txt --recode --recode-INFO-all --out data/uncurated/$c.$f.$t
  mv data/uncurated/$c.$f.$t.recode.vcf data/uncurated/$c.$f.$t.vcf
  grep -v "#" data/uncurated/$c.$f.$t.vcf |wc
done



################################### PCAs #######################################

# PCAs were run on the following sets (and combinations of them):
# Filtered - sites that pass the DUPHOLD and quality filters
# Removed - sites that were removed in the DUPHOLD and quality filtering
# Uncurated - sites in the "Filtered" set but did not have two individuals
#   representing each genotype, and hence could not be manually curated)
# Curated - sites that passed the manual curation in plotcritic
# Rejected - sites that failed the manual curation

mkdir plink/

# PCA in Plink
for type in "DUP" "INV" "DEL"
do
  filt="strict"
  set="curated"
  plink2 --vcf data/$set/1cur.$filt.$type.vcf --dog --make-pgen --sort-vars --out plink/$set.1cur.$filt.$type
  plink2 --pfile plink/$set.1cur.$filt.$type --dog --pca 20 --threads 1 --out plink/$set.1cur.$filt.$type.pca
  plink2 --pfile plink/$set.1cur.$filt.$type --dog --pca 20 --make-rel --threads 1 --out plink/$set.1cur.$filt.$type.pca
done

# Get the diagonal from the rel file, needed for percent explained
for type in "DUP" "INV" "DEL"
do
  awk -v t=$type '{sum+=$NR;}END{print t"\t"sum}' plink/curated.1cur.strict.$type.pca.rel >>plink/final.diagsum
done

# The duphold filtered and uncurated are not curation dependent, so there is
# only one file per type
for type in "DUP" "INV" "DEL"
do
  plink2 --vcf data/filtered/duphold.fold.MHSQ.$type.vcf --dog --make-pgen --sort-vars --out plink/duphold.fold.MHSQ.$type
  plink2 --pfile plink/duphold.fold.MHSQ.$type --dog --pca 20 --threads 1 --out plink/duphold.fold.MHSQ.$type.pca
  plink2 --vcf data/filtered/duphold.removed.$type.vcf --dog --make-pgen --sort-vars --out plink/duphold.removed.$type
  plink2 --pfile plink/duphold.removed.$type --dog --pca 20 --threads 1 --out plink/duphold.removed.$type.pca
  plink2 --vcf data/uncurated/1cur.strict.$type.vcf --dog --make-pgen --sort-vars --out plink/uncurated.$type
  plink2 --pfile plink/uncurated.$type --dog --pca 20 --threads 1 --out plink/uncurated.$type.pca
done

# Combination of rejected and removed (only using 1cur/relaxed)
for type in "DUP" "INV" "DEL"
do
  plink2 --vcf data/rejected_removed/$type.vcf --dog --make-pgen --sort-vars --out plink/rejected_removed.$type
  plink2 --pfile plink/rejected_removed.$type --dog --pca 20 --threads 1 --out plink/rejected_removed.$type.pca
done



############################# SV LENGTH STATS  #################################
# All code for extracting and checking lengths are collected here:
bash/SV_lengths.sh




################################# ADDING SNPs ##################################

# Principal component analysis on SNPs from Smeds and Ellegren 2023.
echo '#!/bin/sh
#plink2 --vcf 100S95F14R.chr1-38.allSNPs.wfm.vcf.gz --dog --make-pgen --sort-vars --out plink/SNPs
plink2 --pfile plink/SNPs --dog --pca 20 --make-rel --threads 4 --out plink/SNPs.pca
' |sbatch -t 1:00:00

# Relation covariance matrix, sum diagonal
awk -v t="SNP" '{sum+=$NR;}END{print t"\t"sum}' plink/SNPs.pca.rel >>plink/final.diagsum




################################### REPEATS ####################################

# Download repeats from UCSC ftp server (annotation from 2011)
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.out.gz
mv canFam3.fa.out.gz reference/repeats.canFam3.fa.out.gz
#Bed file with family names
zcat reference/repeats.canFam3.fa.out.gz |tail -n+4 |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$11}' |sed 's/^chr//' >reference/CanFam3.1.repeats.family.bed
# Bedfile with exact names
zcat reference/repeats.canFam3.fa.out.gz |tail -n+4 |awk -v OFS="\t" '{s=$6-1; print $5,s,$7,$10}' |sed 's/^chr//' >reference/CanFam3.1.repeats.names.bed


# check total length of autosomal repeats:
awk '(/^[1-9]/){sum+=$3-$2}END{print sum}' reference/CanFam3.1.repeats.family.bed
#901692792

# Check if there are overlapping repeats, how many and the sum of the overlap
awk -v OFS="\t" '{if(NR==1){chr=$1; last=$3}else{if(last>$2 && chr==$1){o=last-$2; print line $0,o}; last=$3; chr=$1}}' reference/CanFam3.1.repeats.family.bed |wc
# 95572
awk -v OFS="\t" '{if(NR==1){chr=$1; last=$3}else{if(last>$2 && chr==$1){o=last-$2; print line $0,o}; last=$3; chr=$1}}' reference/CanFam3.1.repeats.family.bed |awk '{sum+=$5}END{print sum}'
# 1349400
# not so many, of the total

# Create bed files of the final SVs
grep -v "#" data/curated/1cur.strict.DEL.vcf |cut -f1-8 |sed 's/SVLEN=-//' |awk -v OFS="\t" '{split($8,t,";"); start=$2-1; e=start+t[2]; print $1,start,e,t[2]}' >data/bed/DEL.final.bed
for type in "DUP" "INV"
do
  grep -v "#" data/curated/1cur.strict.$type.vcf |cut -f1-8 |sed 's/SVLEN=//'|awk -v OFS="\t" '{split($8,t,";"); start=$2-1; e=start+t[2]; print $1,start,e,t[2]}' >data/bed/$type.final.bed
done

# Check overlaps with repeats
intersectBed -wo -a data/bed/DEL.final.bed -b reference/CanFam3.1.repeats.bed |less
# Check number of repeats overlapping with SVs
intersectBed -wo -a data/bed/DEL.final.bed -b reference/CanFam3.1.repeats.bed |wc
# Check how many of each repeats class:
intersectBed -wo -a data/bed/DEL.final.bed -b reference/CanFam3.1.repeats.bed |cut -f7 |sort |uniq -c |less

# Check what repeats are in the peak regions:
# (take an arbirary distance around peaks)
awk '($4>50 && $4<60){print}' data/bed/DEL.final.bed |intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |cut -f7 |sort |uniq -c |less
awk '($4>180 && $4<200){print}' data/bed/DEL.final.bed |intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |cut -f7 |sort |uniq -c |less

# Make file with all repeat families, their lengths and also a class for non-
# repetetive sequence (use custom python script to make sure no bases are
# counted twice, in the case of overlapping repeat annotations)
for type in "DEL" "DUP" "INV"
do
  intersectBed -wo -a data/bed/$type.final.bed -b reference/CanFam3.1.repeats.bed |cut -f1 -d"/" |sed 's/?//' |sed 's/rRNA/RNA/' |sed 's/snRNA/RNA/' |sed 's/tRNA/RNA/' |sed 's/scRNA/RNA/'  >temp.overlap.$type.rep.txt
  python3 python/summarize_repeats.py -f temp.overlap.$type.rep.txt -o data/repeats/$type.repeat.annot.per.base.bed
  cut -f4 data/repeats/$type.repeat.annot.per.base.bed |sort |uniq -c |awk -v OFS="\t" '{print $2,$1}' > data/repeats/$type.repeat.summary.txt
done

# The same as above, but only for the peak regions
# DELETIONS Peak 40-60bp
awk '($4>50 && $4<60){print}' data/bed/DEL.final.bed | intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |cut -f1 -d"/" |sed 's/?//' |sed 's/rRNA/RNA/' |sed 's/snRNA/RNA/' |sed 's/tRNA/RNA/' |sed 's/scRNA/RNA/'  >temp.overlap.DEL.peak50.rep.txt
python3 python/summarize_repeats.py -f temp.overlap.DEL.peak50.rep.txt -o data/repeats/DEL.peak50.repeat.annot.per.base.bed
cut -f4 data/repeats/DEL.peak50.repeat.annot.per.base.bed |sort |uniq -c |awk -v OFS="\t" '{print $2,$1}' > data/repeats/DEL.peak50.repeat.summary.txt
# DELETIONS Peak 180-200bp
awk '($4>180 && $4<200){print}' data/bed/DEL.final.bed | intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |cut -f1 -d"/" |sed 's/?//' |sed 's/rRNA/RNA/' |sed 's/snRNA/RNA/' |sed 's/tRNA/RNA/' |sed 's/scRNA/RNA/'  >temp.overlap.DEL.peak190.rep.txt
python3 python/summarize_repeats.py -f temp.overlap.DEL.peak190.rep.txt -o data/repeats/DEL.peak190.repeat.annot.per.base.bed
cut -f4 data/repeats/DEL.peak190.repeat.annot.per.base.bed |sort |uniq -c |awk -v OFS="\t" '{print $2,$1}' > data/repeats/DEL.peak190.repeat.summary.txt
# DUPLICATIONS Peak 140-180bp
awk '($4>140 && $4<180){print}' data/bed/DUP.final.bed | intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |cut -f1 -d"/" |sed 's/?//' |sed 's/rRNA/RNA/' |sed 's/snRNA/RNA/' |sed 's/tRNA/RNA/' |sed 's/scRNA/RNA/'  >temp.overlap.DUP.peak160.rep.txt
python3 python/summarize_repeats.py -f temp.overlap.DUP.peak160.rep.txt -o data/repeats/DUP.peak160.repeat.annot.per.base.bed
cut -f4 data/repeats/DUP.peak160.repeat.annot.per.base.bed |sort |uniq -c |awk -v OFS="\t" '{print $2,$1}' > data/repeats/DUP.peak160.repeat.summary.txt


# For large Deletion peak;
#check peak variants with multiple annotations and check how many LINE,
awk '($4>180 && $4<200){print}' data/bed/DEL.final.bed | intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |cut -f2 |uniq -c |awk '($1>1){print $2}' |grep -f - temp.overlap.DEL.peak190.rep.txt |grep "LINE" |wc
#2078
# Then check number of SINE overlapping variants in total.
awk '($4>180 && $4<200){print}' data/bed/DEL.final.bed | intersectBed -wo -a - -b reference/CanFam3.1.repeats.bed |grep "SINE" |cut -f2 |uniq |wc
# 10176



########################## OUTGROUPS FOR POLARIZATION ##########################
# 2023-09-24

# Download outgroups
job=1000000
for ind in $(grep -v "#" helpfiles/SRA_accession_outgroups.txt |cut -f2)
do
  echo $ind
  sbatch -J $ind.wget -p core  -t 10:00:00  bash/run_ENA_download.sh $ind helpfiles/SRA_accession_outgroups.txt |cut -f4 -d" ")
done

# Map reads with snakemake pipeline (code from Smeds & Ellegren et al 2023)
snakemake --snakefile snakemake/outgroup.snakefile -p -j 64  --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -e {cluster.error} -o {cluster.output} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}" --cluster-config snakemake/cluster_outgroup.json

# Genotype variants in outgroups (this was run with Lars smoove code, and
# therefore a later version of Snakemake: snakemake/7.18.2)
mkdir -p smoove/genotyped
# run SV caller in snakemake
snakemake -np --snakefile snakemake/smoove.snakefile
snakemake -s snakemake/smoove.snakefile --use-singularity --nolock -j 64 \
        --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} -J {cluster.name} --parsable -e {cluster.error}.slurm -o {cluster.output}.slurm --mail-user {cluster.email} --mail-type {cluster.mailtype}" \
        --cluster-config snakemake/cluster_smoove.json

# Extract relevant columns
for i in "AfGWo" "BlBJa"
do
  zcat smoove/genotyped/$i.joint-smoove.genotyped.vcf.gz |grep -v "##" |cut -f1,2,3,10  |cut -f1 -d":" >smoove/genotyped/$i.extract.txt
done

# Polarize all variants
for i in "AfGWo" "BlBJa"
do
  grep -v "_" smoove/genotyped/$i.extract.txt |cut -f3,4 |sort >smoove/genotyped/$i.all.txt
done
# Then join the two outgroups and print the posititons where they agree
paste smoove/genotyped/AfGWo.all.txt smoove/genotyped/BlBJa.all.txt |sed 's/ /\t/g' |awk '($1==$3){if($2==$4 && $2!="0/1"){print}}' |grep -v '\./\.' |cut -f1,2  >data/polarization/joint.all.txt

# Check length stats for each class
for type in "DEL" "DUP" "INV"
do
  echo $type
  join <(sort data/polarization/$type.joint.final.txt) <(less data/lengths/curated_lengths_with_ID.txt |awk '{split($1,s,";"); print s[1]"\t"$2}' |sort) >$type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3<=50){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>50 && $3<=100){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>100 && $3<=200){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>200 && $3<=300){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>300 && $3<=400){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>400 && $3<=500){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>500 && $3<=1000){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>1000 && $3<=10000){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>10000 && $3<=50000){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
  awk -v a0=0 -v a1=0 '($3>50000){n++; if($2=="0/0"){a0=a0+1}else if($2=="1/1"){a1=a1+1}}END{print "Ancestral 0/0: "a0" Ancestral 1/1: "a1" Total: "n}' $type.tmp.pol.txt
done


############################## VEP ON FINAL SITES ##############################
# NOTE: The VEP impact output (HIGH, LOW, MODERATE, MODIFIER) was not used in
# the final study, but we extracted overlap with genes from its output.

# Run VEP (code taken from Smeds & Ellegren 2023)
for type in "DUP" "INV" "DEL"
do
  echo $type
echo '#!/bin/sh
vep --cache --dir $VEP_CACHE -i data/curated/1cur.strict.'$type'.vcf -o vep/'$type'.raw_output.txt --species canis_familiaris --fork 4 --force_overwrite --sift b
' | sbatch -t 2-00:00:00 -n4 -J $type.vep
done

# Also for duphold filtered deletions 50-10000bp
echo '#!/bin/sh
vep --cache --dir $VEP_CACHE -i data/filtered/duphold.fold.MHSQ.DEL.50-10000.vcf -o vep/DEL.duphold50-10000.raw_output.txt --species canis_familiaris --fork 4 --force_overwrite --sift b
' | sbatch  -t 5:00:00 -n4 -J duphold.vep

# Code for extracting SVs in coding, intronic and intergenic regions and count
# number of variants for tables are found in:
bash/genes.sh



################################# TRIO CHECK ###################################
# Checking concordance between parents and offspring in 7 different trios

# Making vcf files with the trios
for type in "DEL" "DUP" "INV"
do
    # create trio subsets
    bcftools view -s 83-G31-13,82-G23-13,104-G100-14 data/curated/1cur.strict.${type}.vcf > data/trios/curated.se1.${type}.vcf
    bcftools view -s 83-G31-13,82-G23-13,89-G67-15 data/curated/1cur.strict.${type}.vcf > data/trios/curated.se2.${type}.vcf
    bcftools view -s 6-M-98-02,11-M-98-03,13-M-98-08 data/curated/1cur.strict.${type}.vcf > data/trios/curated.se3.${type}.vcf
    bcftools view -s 6-M-98-02,11-M-98-03,Varg0109 data/curated/1cur.strict.${type}.vcf > data/trios/curated.se4.${type}.vcf
    bcftools view -s W57r,W45,W58 data/curated/1cur.strict.${type}.vcf > data/trios/curated.fi1.${type}.vcf
    bcftools view -s W53r,W54r,W60 data/curated/1cur.strict.${type}.vcf > data/trios/curated.fi2.${type}.vcf
    bcftools view -s W53r,W54r,W62r data/curated/1cur.strict.${type}.vcf > data/trios/curated.fi3.${type}.vcf
done

# count how many of the sites have matching genotypes (or not)
for type in "DEL" "DUP" "INV"
do
  for trio in "se1" "se2" "se3" "se4" "fi1" "fi2" "fi3"
  do
    bcftools view -h data/trios/curated.$trio.$type.vcf | tail -n1 | cut -f 10- |awk '{print $0"\tNumber"}' > data/trios/$trio.$type.GTcount.txt
    grep -v "#" data/trios/curated.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' | sort |uniq -c |awk -v OFS="\t" '{print $2,$3,$4,$1}' >>data/trios/$trio.$type.GTcount.txt
    bad=`grep -v "#" data/trios/curated.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -e $'0/0 0/0 0/1' -e $'0/0 0/0 1/1' -e $'0/0 0/1 1/1' -e $'0/0 1/1 0/0' -e $'0/0 1/1 1/1' -e $'0/1 0/0 1/1' -e $'0/1 1/1 0/0' -e $'1/1 0/0 0/0' -e $'1/1 0/0 1/1' -e $'1/1 0/1 0/0' -e $'1/1 1/1 0/0' -e $'1/1 1/1 0/1' |wc -l`
    ok=`grep -v "#" data/trios/curated.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -v "\./\." |grep -v -e $'0/0 0/0 0/1' -e $'0/0 0/0 1/1' -e $'0/0 0/1 1/1' -e $'0/0 1/1 0/0' -e $'0/0 1/1 1/1' -e $'0/1 0/0 1/1' -e $'0/1 1/1 0/0' -e $'1/1 0/0 0/0' -e $'1/1 0/0 1/1' -e $'1/1 0/1 0/0' -e $'1/1 1/1 0/0' -e $'1/1 1/1 0/1' |wc -l`
    echo "$type $trio $ok $bad" |sed 's/ /\t/g'
  done
done

# Some more code for making example figures, and looking at concordance in the
# SNP data from Smeds & Ellegren, 2023.
bash/trios.sh



############################ PLOTTING AND STATS IN R ###########################

# Plot with the length distributions for the different classes
R plot_figure1.R
# Add uncurated and rejected/removed
R plot_figureS3.R

# Plot repeat content in SVs
R plot_figure2.R

# Plot PCAs of SVs and SNPs
R plot_figure3.R
# PCA comparing one and two curators
R plot_figureS1.R
# PCA comparing kept with removed and uncurated
R plot_figureS2.R

# Plot Allele frequency spectrum (using deletions before manual curation, to
# not lose rare variants)
R plot_figure4.R

# Plot load for each generation
R plot_figure5.R

# Plot load for descendants to original founders compared to descendants to
# immigrants from the same time period
R plot_figure6.R

