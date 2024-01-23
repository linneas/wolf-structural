# COMMANDS FOR STRUCTURAL VARIATION IN FENNOSCANDIAN VOLVES
# written by Linnéa Smeds and Lars S.A. Huson


# Software versions used:
#vcftools/0.1.16
#bcftools/1.17
#snakemake/7.18.2

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
# deletions (DHFFC<0.7):
bcftools view -i 'SVTYPE="DEL" & DHFFC<0.7' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf >data/filtered/duphold.fold.DEL.vcf
# 61476 DEL
# duplications (DHFFC>1.3):
bcftools view -i 'SVTYPE="DUP" & DHFFC>1.3' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' data/smoove/annotated/cohort.anno.filtered.GT.vcf >data/filtered/duphold.fold.DUP.vcf
# 6645 DUP
# inversions (no filter):
cp data/smoove/annotated/cohort.anno.filtered.GT.vcf data/filtered/duphold.fold.INV.vcf
# 9077 INV
