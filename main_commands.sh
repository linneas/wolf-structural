# COMMANDS FOR STRUCTURAL VARIATION IN FENNOSCANDIAN VOLVES
# written by Linnéa Smeds and Lars S.A. Huson


#### Software versions used:
#vcftools/0.1.16
#bcftools/1.17
#snakemake/7.18.2
#python/3.9.5
#pysam/0.17.0-python3.9.5

############################## SMOVE AND DUPHOLD ###############################
#Code for running smoove:
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
# Save only variants >50bp and <10,000bp
bcftools view -i 'SVLEN>=-10000 & SVLEN<-50' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/filtered/duphold.fold.MHSQ.DEL.vcf >data/filtered/duphold.fold.MHSQ.DEL.50-10000.vcf



########################### SAMPLOT AND PLOTCRITIC #############################
# (Deletions are used as example)
# The variants are divided by length to simply the manual curation

# Preparation, create multiple scripts for each batch to speed up runtime
bash/prepare_samplot_scripts_del.sh

# Samplot was installed in a conda environment
module unload python pysam
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
### genertate plotcritic websites (one for each length batch)
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
  # Get the rejected sites and ids (first strict, then relaxed)
  grep -e $'0.0\t100.0\t0.0' -e $'0.0\t0.0\t100.0' plotcritic/plcr_${t}_summary_report_Linnea.tsv |awk '{s=$3+1; print $2"\t"s"\t"$4}' >plotcritic/rejected_${t}_1cur_$f.pos.txt
  awk '{print $1":"$2"-"$3}' plotcritic/rejected_${t}_1cur_$f.pos.txt |sort |join - <(sort data/filtered/duphold.fold.MHSQ.$t.list) |cut -f2 -d" " |sort -g >plotcritic/rejected_${t}_1cur_$f.id.txt
  vcftools --vcf data/filtered/duphold.fold.MHSQ.$t.vcf --snps plotcritic/rejected_${t}_1cur_${f}.id.txt --recode --recode-INFO-all --out data/rejected/1cur.$f.$t
  mv data/rejected/1cur.$f.$t.recode.vcf data/rejected/1cur.$f.$t.vcf
done
