#!/bin/bash

# Code by Lars Huson and Linnéa Smeds

# Divided the deletions in batches based on length
mkdir data/filtered/del_by_lengths/
bcftools view -i 'SVLEN>-50 & SVLEN<=0' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.0-50.vcf
bcftools view -i 'SVLEN>=-100 & SVLEN<-50' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.50-100.vcf
bcftools view -i 'SVLEN>=-150 & SVLEN<-100' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.100-150.vcf
bcftools view -i 'SVLEN>=-175 & SVLEN<-150' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.150-175.vcf
bcftools view -i 'SVLEN>=-200 & SVLEN<-175' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.175-200.vcf
bcftools view -i 'SVLEN>=-250 & SVLEN<-200' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.200-250.vcf
bcftools view -i 'SVLEN>=-300 & SVLEN<-250' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.250-300.vcf
bcftools view -i 'SVLEN>=-400 & SVLEN<-300' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.300-400.vcf
bcftools view -i 'SVLEN>=-500 & SVLEN<-400' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.400-500.vcf
bcftools view -i 'SVLEN>=-750 & SVLEN<-500' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.500-750.vcf
bcftools view -i 'SVLEN>=-1000 & SVLEN<-750' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.750-1000.vcf
bcftools view -i 'SVLEN<=-1000' -f '%CHROM\t%POS\t%INFO/END[\t%GT]\n' smoove/annotated/cohort.anno.filtered.del.vcf >data/filtered/del_by_length/cohort.anno.filtered.del.from1000.vcf

#write bash code with samplot commands
# The script gen_samplot_del.py was kindly proivided by Gabriel David 
for set in "0-50" "50-100"  "100-150" "150-175" "175-200" "200-250" "250-300" "300-400" "400-500" "500-750" "750-1000" "from1000"  #
do
  mkdir samplot/imgs_DEL_$set
  python python/gen_samplot_del.py --vcf data/filtered/del_by_length/cohort.anno.filtered.del.$set.vcf --contigs helpfiles/contigs.json --bams helpfiles/bamlocations.json --outdir samplot/imgs_DEL_$set >batchscripts/samplot_dup_cmds.$set.sh
done

#To run it faster on a slurm cluster, split into multiple slurm files
mkdir batchscripts/parts
mkdir batchscripts/tmp
for set in "0-50" "50-100" "100-150" "150-175" "175-200" "200-250" "250-300" "300-400" "400-500" "500-750" "750-1000" "from1000"
do
  split -d -l 100 batchscripts/samplot_dup_cmds.$set.sh batchscripts/tmp/samplot_dup_cmds.$set.

  for file in $(ls batchscripts/tmp/samplot_dup_cmds.$set.*)
  do
    ls $file
    i=`echo $file |cut -f3 -d"."`
    out=`echo $file |sed 's/tmp/parts/' `
    echo $out
    echo '#!/bin/bash
#SBATCH -p core -n 1
#SBATCH -t 1:00:00
#SBATCH -J samplot_del_'$set'_'$i'
' >$out.sh
      cat $file >>$out.sh
  done
done
