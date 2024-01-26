# Code by LS and LSAH


# Some illustrative examples for the paper;
# one with homozygous parents heterozygous offsp.,
# and one with heterozygous parents and homozygous offspring.
# Use Scandinavian family
for type in "DEL" "DUP" "INV"
do
  # create family
  bcftools view -s 6-M-98-02,11-M-98-03,13-M-98-08,Varg0109 data/curated/1cur.strict.${type}.vcf > data/trios/curated.family2.${type}.vcf
  grep -v "#" data/trios/curated.family2.$type.vcf |cut -f 1-3,10- |awk '{for(i=4; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' >data/trios/curated.family2.${type}.GT.txt
done
# Parents homozygous (take a deletion)
awk '((($4=="0/0" && $5=="1/1") || ($4=="1/1" && $5=="0/0")) && ($6=="0/1" && $7=="0/1")){print}' data/trios/curated.family2.DEL.GT.txt |join -1 3 -2 1 - <(sort data/lengths/curated_lengths_with_ID.txt)
# Parents heterozygous (take a duplication)
awk '(($4=="0/1" && $5=="0/1") && (($6=="1/1" && $7=="0/0") || ($6=="0/0" && $7=="1/1"))){print}' data/trios/curated.family2.DUP.GT.txt |join -1 3 -2 1 - <(sort data/lengths/curated_lengths_with_ID.txt)

# Plot in the same way as before
conda activate samplotcritic
samplot plot -r data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa -n 6-M-98-02 11-M-98-03 13-M-98-08 Varg0109 -b data/raw_data/cram/6-M-98-02.cram data/raw_data/cram/11-M-98-03.cram data/raw_data/cram/13-M-98-08.cram data/raw_data/cram/Varg0109.cram -o samplot/family_test/DEL_1_87552894_87554098_171_192_90_55_202_105.png -c 1 -s 87552894 -e 87554098 -a  -t DEL
samplot plot -r data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa -n Parent1 Parent2 Offspring1 Offspring2 -b data/raw_data/cram/6-M-98-02.cram data/raw_data/cram/11-M-98-03.cram data/raw_data/cram/13-M-98-08.cram data/raw_data/cram/Varg0109.cram -o samplot/family_test/DEL_1_87552894_87554098_171_192_90_55_202_105.png -c 1 -s 87552894 -e 87554098 -a  -t DEL
samplot plot -r data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa -n Parent1 Parent2 Offspring -b data/raw_data/cram/6-M-98-02.cram data/raw_data/cram/11-M-98-03.cram data/raw_data/cram/13-M-98-08.cram -o samplot/family_test/DEL_1_87552894_87554098_1Offpr.png -c 1 -s 87552894 -e 87554098 -a  -t DEL
samplot plot -r data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa -n Parent1 Parent2 Offspring1 Offspring2 -b data/raw_data/cram/6-M-98-02.cram data/raw_data/cram/11-M-98-03.cram data/raw_data/cram/13-M-98-08.cram data/raw_data/cram/Varg0109.cram -o samplot/family_test/DUP_8_32844260_32844773_192_59_135_173_159_111.png -c 8 -s 32844260 -e 32844773 -a  -t DUP

# COMPARE CONCORDANCE FOR SNPs
# create trio subsets
bcftools view -s 83-G31-13,82-G23-13,104-G100-14 /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/se1.SNPs.vcf
bcftools view -s 83-G31-13,82-G23-13,89-G67-15 /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/se2.SNPs.vcf
bcftools view -s 6-M-98-02,11-M-98-03,13-M-98-08 /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/se3.SNPs.vcf
bcftools view -s 6-M-98-02,11-M-98-03,Varg0109 /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/se4.SNPs.vcf
bcftools view -s W57r,W45,W58 /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/fi1.SNPs.vcf
bcftools view -s W53r,W54r,W60 /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/fi2.SNPs.vcf
bcftools view -s W53r,W54r,W62r /proj/sllstore2017034/repos/Deleterious_mutations/vcf/Pol.2out.mac2/100S95F14R.chr1-38.allSNPs.wfm.vcf.gz > data/trios/fi3.SNPs.vcf

for trio in "fi2" "fi3" "se1" "se2" "se3" "fi1" "fi2" "fi3" "se4"
do
    bad=`grep -v "#" data/trios/$trio.SNPs.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -e $'0/0 0/0 0/1' -e $'0/0 0/0 1/1' -e $'0/0 0/1 1/1' -e $'0/0 1/1 0/0' -e $'0/0 1/1 1/1' -e $'0/1 0/0 1/1' -e $'0/1 1/1 0/0' -e $'1/1 0/0 0/0' -e $'1/1 0/0 1/1' -e $'1/1 0/1 0/0' -e $'1/1 1/1 0/0' -e $'1/1 1/1 0/1' |wc -l`
    ok=`grep -v "#" data/trios/$trio.SNPs.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -v "\./\." |grep -v -e $'0/0 0/0 0/1' -e $'0/0 0/0 1/1' -e $'0/0 0/1 1/1' -e $'0/0 1/1 0/0' -e $'0/0 1/1 1/1' -e $'0/1 0/0 1/1' -e $'0/1 1/1 0/0' -e $'1/1 0/0 0/0' -e $'1/1 0/0 1/1' -e $'1/1 0/1 0/0' -e $'1/1 1/1 0/0' -e $'1/1 1/1 0/1' |wc -l`
    echo "$trio $ok $bad" |sed 's/ /\t/g'
  done
done


# COMPARE WITH HOW IT LOOKS BEFORE CURATION
for type in "DEL" "DUP" "INV"
do
    # create trio subsets
    bcftools view -s 83-G31-13,82-G23-13,104-G100-14 data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.se1.${type}.vcf
    bcftools view -s 83-G31-13,82-G23-13,89-G67-15 data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.se2.${type}.vcf
    bcftools view -s 6-M-98-02,11-M-98-03,13-M-98-08 data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.se3.${type}.vcf
    bcftools view -s 6-M-98-02,11-M-98-03,Varg0109 data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.se4.${type}.vcf
    bcftools view -s W57r,W45,W58 data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.fi1.${type}.vcf
    bcftools view -s W53r,W54r,W60 data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.fi2.${type}.vcf
    bcftools view -s W53r,W54r,W62r data/filtered/duphold.fold.MHSQ.${type}.vcf > data/trios/dupfilt.fi3.${type}.vcf
done

for type in "DEL" "DUP" "INV"
do
  for trio in "se1" "se2" "se3" "se4" "fi1" "fi2" "fi3"
  do
    bcftools view -h data/trios/curated.$trio.$type.vcf | tail -n1 | cut -f 10- |awk '{print $0"\tNumber"}' > data/trios/$trio.$type.GTcount.txt
    grep -v "#" data/trios/curated.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' | sort |uniq -c |awk -v OFS="\t" '{print $2,$3,$4,$1}' >>data/trios/$trio.$type.GTcount.txt
    bad=`grep -v "#" data/trios/dupfilt.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -e $'0/0 0/0 0/1' -e $'0/0 0/0 1/1' -e $'0/0 0/1 1/1' -e $'0/0 1/1 0/0' -e $'0/0 1/1 1/1' -e $'0/1 0/0 1/1' -e $'0/1 1/1 0/0' -e $'1/1 0/0 0/0' -e $'1/1 0/0 1/1' -e $'1/1 0/1 0/0' -e $'1/1 1/1 0/0' -e $'1/1 1/1 0/1' |wc -l`
    ok=`grep -v "#" data/trios/dupfilt.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -v "\./\." |grep -v -e $'0/0 0/0 0/1' -e $'0/0 0/0 1/1' -e $'0/0 0/1 1/1' -e $'0/0 1/1 0/0' -e $'0/0 1/1 1/1' -e $'0/1 0/0 1/1' -e $'0/1 1/1 0/0' -e $'1/1 0/0 0/0' -e $'1/1 0/0 1/1' -e $'1/1 0/1 0/0' -e $'1/1 1/1 0/0' -e $'1/1 1/1 0/1' |wc -l`
    echo "$type $trio $ok $bad" |sed 's/ /\t/g'
  done
done

# Save to file 
for type in "DEL" "DUP" "INV"
do
  for trio in "se1" "se2" "se3" "se4" "fi1" "fi2" "fi3"
  do
    grep -v "#" data/trios/dupfilt.$trio.$type.vcf |cut -f 10- |awk '{for(i=1; i<=NF; i++){split($i,s,":"); $i=s[1]}; print}' |grep -v "\./\." |awk '{if($0=="0/0 0/0 0/1" || $0=="0/0 0/0 1/1" ||  $0=="0/0 0/1 1/1" || $0=="0/0 1/1 0/0" || $0=="0/0 1/1 1/1" || $0=="0/1 0/0 1/1" || $0=="0/1 1/1 0/0" || $0=="1/1 0/0 0/0" || $0=="1/1 0/0 1/1" || $0=="1/1 0/1 0/0" || $0=="1/1 1/1 0/0" || $0=="1/1 1/1 0/1"){print "0"}else{print "1"}}' >data/trios/dupfilt.$trio.$type.txt
  done
done
