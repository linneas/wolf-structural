#!/bin/bash
# Code by LS



# Make tables with SV types and lengths
grep -v "#" data/curated/1cur.strict.DEL.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=-//' |awk '{print "Deletions\tCurated\t"$1}' >data/lengths/all_lengths.txt
grep -v "#" data/curated/1cur.strict.DUP.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Duplications\tCurated\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/curated/1cur.strict.INV.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Inversions\tCurated\t"$1}' >>data/lengths/all_lengths.txt

grep -v "#" data/uncurated/1cur.strict.DEL.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=-//' |awk '{print "Deletions\tUncurated\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/uncurated/1cur.strict.DUP.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Duplications\tUncurated\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/uncurated/1cur.strict.INV.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Inversions\tUncurated\t"$1}' >>data/lengths/all_lengths.txt

grep -v "#" data/rejected/1cur.strict.DEL.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=-//' |awk '{print "Deletions\tRejected\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/rejected/1cur.strict.DUP.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Duplications\tRejected\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/rejected/1cur.strict.INV.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Inversions\tRejected\t"$1}' >>data/lengths/all_lengths.txt

grep -v "#" data/filtered/duphold.removed.DEL.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=-//' |awk '{print "Deletions\tRemoved\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/filtered/duphold.removed.DUP.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Duplications\tRemoved\t"$1}' >>data/lengths/all_lengths.txt
grep -v "#" data/filtered/duphold.removed.INV.vcf |cut -f8 |cut -f2 -d";"|sed 's/SVLEN=//' |awk '{print "Inversions\tRemoved\t"$1}' >>data/lengths/all_lengths.txt

# Find min and max lengths 
for type in "Deletions" "Duplications" "Inversions"
do
  for set in "Uncurated" "Curated" "Rejected" "Removed"
  do
    echo $type" "$set
    grep $type data/lengths/all_lengths.txt |grep $set | awk -v max=0 -v min=10000000 '{if($3>max){max=$3}; if($3<min){min=$3}}END{print "Min: "min" Max:"max}'
  done
done

# Get length classes for table
for type in "DEL" "DUP" "INV"
do
  echo $type
  for set in "raw_calls" "duphold.fold.MHSQ" "genotype_filtered" "curated"
  do
    echo $set
    awk '($1<=50){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>50 && $1<=100){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>100 && $1<=200){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>200 && $1<=300){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>300 && $1<=400){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>400 && $1<=500){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>500 && $1<=1000){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>1000 && $1<=10000){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>10000 && $1<=50000){print}' data/lengths/$set.$type.lengths |wc -l
    awk '($1>50000){print}' data/lengths/$set.$type.lengths |wc -l
  done
done

# Make a table with lengths and ID
rm -f data/lengths/curated_lengths_with_ID.txt
for type in "DEL" "DUP" "INV"
do
  grep -v "#" data/curated/1cur.strict.$type.vcf |cut -f3,8 |awk '{split($2,s,";"); print $1"\t"s[2]}'|sed 's/SVLEN=//' |sed 's/-//' >>data/lengths/curated_lengths_with_ID.txt
done
