# Code by LS

# Extract genes overlapping with SVs
for type in "DEL" "DUP" "INV"
do
  #All genes affected by SVs
  grep -v "#" vep/$type.raw_output.txt |cut -f4,7,14 |cut -f1 -d";" |grep -v "intergenic_variant" |grep -v "upstream_gene_variant" |grep -v "downstream_gene_variant" |cut -f1 |sort |uniq >vep/$type.all_affected.genes.txt
  # Genes affected by inframe variants
  grep -v "#" vep/$type.raw_output.txt |cut -f4,7,14 |cut -f1 -d";" |grep "inframe" |cut -f1 |sort |uniq >vep/$type.inframe_affected.genes.txt
  # Genes affected by frameshifts variants
  grep -v "#" vep/$type.raw_output.txt |cut -f4,7,14 |cut -f1 -d";" |grep "frameshift" |cut -f1 |sort |uniq >vep/$type.frameshift_affected.genes.txt
  # Genes affected by start/stop overlap
  grep -v "#" vep/$type.raw_output.txt |cut -f4,7,14 |cut -f1 -d";" |grep "stop_lost\|start_lost\|stop_gained" |cut -f1 |sort |uniq >vep/$type.startstop_affected.genes.txt
  # Full genes overlap
  grep -v "#" vep/$type.raw_output.txt |cut -f4,7,14 |cut -f1 -d";" |grep "transcript_ablation\|transcript_amplification" |cut -f1 |sort |uniq >vep/$type.fulloverlap_affected.genes.txt
  # All other CDS overlap (not start/stop, inframe)
  grep -v "#" vep/$type.raw_output.txt |cut -f4,7,14 |cut -f1 -d";" |grep "coding_sequence" |grep -v "start_lost" |grep -v "stop_lost" |cut -f1 |sort |uniq >vep/$type.CDSother_affected.genes.txt
done
# Why is "coding_sequence_variant" only classified as MODIFIER????

# Most severe, in descending order:
# Full overlap, Start/stop, Frameshift, CDSother, Inframe
# When counting genes, I remove overlap with more severe categories
for type in "DEL" "DUP" "INV"
do
  echo $type
  stst=`join -v 1 vep/$type.startstop_affected.genes.txt vep/$type.fulloverlap_affected.genes.txt |wc -l`
  echo "start/stop affected "$stst
  frm=`join -v 1 vep/$type.frameshift_affected.genes.txt vep/$type.fulloverlap_affected.genes.txt |join -v 1 - vep/$type.startstop_affected.genes.txt  |wc -l`
  echo "frameshift affected "$frm
  inf=`join -v 1 vep/$type.inframe_affected.genes.txt vep/$type.fulloverlap_affected.genes.txt |join -v 1 - vep/$type.startstop_affected.genes.txt |join -v 1 - vep/$type.frameshift_affected.genes.txt  |wc -l`
  echo "inframe affected "$inf
  other=`join -v 1 vep/$type.CDSother_affected.genes.txt vep/$type.fulloverlap_affected.genes.txt |join -v 1 - vep/$type.startstop_affected.genes.txt |join -v 1 - vep/$type.frameshift_affected.genes.txt  |join -v 1 - vep/$type.inframe_affected.genes.txt |wc -l`
  echo "Other CDS affected "$other
  nonc=`join -v 1 vep/$type.all_affected.genes.txt vep/$type.fulloverlap_affected.genes.txt |join -v 1 - vep/$type.startstop_affected.genes.txt |join -v 1 - vep/$type.frameshift_affected.genes.txt  |join -v 1 - vep/$type.inframe_affected.genes.txt |join -v 1 - vep/$type.CDSother_affected.genes.txt |wc -l`
    echo "Only non-coding overlap "$nonc
done

# And total number (if there are overlaps between different categories)
full=`cat vep/*.fulloverlap_affected.genes.txt |wc -l`
echo "full overlap affected "$full
stst=`join -v 1 <(cat vep/*.startstop_affected.genes.txt |sort) <(cat vep/*.fulloverlap_affected.genes.txt |sort) |wc -l`
echo "start/stop affected "$stst
frm=`join -v 1 <(cat vep/*.frameshift_affected.genes.txt |sort) <(cat vep/*.fulloverlap_affected.genes.txt|sort) |join -v 1 - <(cat vep/*.startstop_affected.genes.txt |sort)  |wc -l`
echo "frameshift affected "$frm
inf=`join -v 1 <(cat vep/*.inframe_affected.genes.txt |sort) <(cat vep/*.fulloverlap_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.startstop_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.frameshift_affected.genes.txt |sort)  |wc -l`
echo "inframe affected "$inf
other=`join -v 1 <(cat vep/*.CDSother_affected.genes.txt |sort) <(cat vep/*.fulloverlap_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.startstop_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.frameshift_affected.genes.txt |sort)  |join -v 1 - <(cat vep/*.inframe_affected.genes.txt |sort) |wc -l`
echo "Other CDS affected "$other
nonc=`join -v 1 <(cat vep/*.all_affected.genes.txt |sort) <(cat vep/*.fulloverlap_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.startstop_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.frameshift_affected.genes.txt |sort)  |join -v 1 - <(cat vep/*.inframe_affected.genes.txt |sort) |join -v 1 - <(cat vep/*.CDSother_affected.genes.txt |sort) |wc -l`
  echo "Only non-coding overlap "$nonc

# Extract variants instead of genes
for type in "DEL" "DUP" "INV"
do
 grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "transcript_ablation\|transcript_amplification\|start_lost\|stop_lost\|frameshift\|inframe\|coding_sequence_variant" |cut -f1 |sort |uniq  >vep/$type.deleterious.txt
done

# Get a list of intergenic and intronic variants, that are NOT also in the genic
# set (= make such a set first)
for type in "INV" #"DEL" "DUP" "INV"
do
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep -v "intergenic" |grep "transcript_ablation\|transcript_amplification\|start_lost\|stop_lost\|frameshift\|inframe\|coding" |cut -f1 |sort |uniq >vep/$type.some_gene_overlap.txt
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "intron" |cut -f1 |sort |uniq |join -v 1 - vep/$type.some_gene_overlap.txt >vep/$type.only_intron.txt
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "intergenic" |cut -f1 |sort |uniq |join -v 1 - vep/$type.some_gene_overlap.txt |join -v1 - vep/$type.only_intron.txt >vep/$type.only_intergenic.txt
done

# Get number of variants for all categories, for table
for type in "DEL" "DUP" "INV"
do
  # All SVs affecting genes
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep -v "intergenic_variant" |grep -v "upstream_gene_variant" |grep -v "downstream_gene_variant" |cut -f1 |sort |uniq >vep/$type.affecting_genes.txt
  # Genes affected by inframe variants
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "inframe" |cut -f1 |sort |uniq >vep/$type.inframe_affecting.txt
  # Genes affected by frameshifts variants
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "frameshift" |cut -f1 |sort |uniq >vep/$type.frameshift_affecting.txt
  # Genes affected by start/stop overlap
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "stop_lost\|start_lost\|stop_gained" |cut -f1 |sort |uniq >vep/$type.startstop_affecting.txt
  # Full genes overlap
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "transcript_ablation\|transcript_amplification" |cut -f1 |sort |uniq >vep/$type.fulloverlap_affecting.txt
  # All other CDS overlap (not start/stop, inframe)
  grep -v "#" vep/$type.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "coding_sequence" |grep -v "start_lost" |grep -v "stop_lost" |cut -f1 |sort |uniq >vep/$type.CDSother_affecting.txt
done

# More numbers:
for type in "DEL" "DUP" "INV"
do
  echo $type
  full=`cat vep/$type.fulloverlap_affecting.txt |wc -l`
  echo "full gene overlap affecting "$full
  stst=`join -v 1 vep/$type.startstop_affecting.txt vep/$type.fulloverlap_affecting.txt |wc -l`
  echo "start/stop affecting "$stst
  frm=`join -v 1 vep/$type.frameshift_affecting.txt vep/$type.fulloverlap_affecting.txt |join -v 1 - vep/$type.startstop_affecting.txt  |wc -l`
  echo "frameshift affecting "$frm
  inf=`join -v 1 vep/$type.inframe_affecting.txt vep/$type.fulloverlap_affecting.txt |join -v 1 - vep/$type.startstop_affecting.txt |join -v 1 - vep/$type.frameshift_affecting.txt  |wc -l`
  echo "inframe affecting "$inf
  other=`join -v 1 vep/$type.CDSother_affecting.txt vep/$type.fulloverlap_affecting.txt |join -v 1 - vep/$type.startstop_affecting.txt |join -v 1 - vep/$type.frameshift_affecting.txt  |join -v 1 - vep/$type.inframe_affecting.txt |wc -l`
  echo "Other CDS affecting "$other
  nonc=`join -v 1 vep/$type.affecting_genes.txt vep/$type.fulloverlap_affecting.txt |join -v 1 - vep/$type.startstop_affecting.txt |join -v 1 - vep/$type.frameshift_affecting.txt  |join -v 1 - vep/$type.inframe_affecting.txt |join -v 1 - vep/$type.CDSother_affecting.txt |wc -l`
    echo "Only non-coding overlap "$nonc
done

# Get some stats for the uncurated set of deletions 50-10,000
grep -v "#" vep/DEL.duphold50-10000.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "transcript_ablation\|transcript_amplification\|start_lost\|stop_lost\|frameshift\|inframe\|coding_sequence_variant" |cut -f1 |sort |uniq  >vep/DEL.duphold50-10000.deleterious.txt
grep -v "#" vep/DEL.duphold50-10000.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep -v "intergenic" |grep "transcript_ablation\|transcript_amplification\|start_lost\|stop_lost\|frameshift\|inframe\|coding" |cut -f1 |sort |uniq >vep/DEL.duphold50-10000.some_gene_overlap.txt
grep -v "#" vep/DEL.duphold50-10000.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "intron" |cut -f1 |sort |uniq |join -v 1 - vep/DEL.duphold50-10000.some_gene_overlap.txt >vep/DEL.duphold50-10000.only_intron.txt
grep -v "#" vep/DEL.duphold50-10000.raw_output.txt |cut -f1,4,7,14 |cut -f1 -d";" |grep "intergenic" |cut -f1 |sort |uniq |join -v 1 - vep/DEL.duphold50-10000.some_gene_overlap.txt |join -v1 - vep/DEL.duphold50-10000.only_intron.txt >vep/DEL.duphold50-10000.only_intergenic.txt

#Merge these to one single vep file
for type in "DEL" "DUP" "INV" "DEL.duphold50-10000"
do
  awk '{print $1"\tdeleterious"}' vep/$type.deleterious.txt |cat - <(awk '{print $1"\tintronic"}' vep/$type.only_intron.txt) |cat - <(awk '{print $1"\tintergenic"}' vep/$type.only_intergenic.txt) >vep/$type.comb.txt
done

# Merge more specific categories (for table)
for type in "DEL" "DUP" "INV"
do
  echo $type
  cat vep/$type.startstop_affecting.txt vep/$type.frameshift_affecting.txt vep/$type.inframe_affecting.txt vep/$type.CDSother_affecting.txt |sort |uniq | join -v 1 - vep/$type.fulloverlap_affecting.txt >vep/$type.partialCDS_affecting.txt
  join -v 1 vep/$type.affecting_genes.txt vep/$type.fulloverlap_affecting.txt |join -v 1 - vep/$type.partialCDS_affecting.txt |join -v 1 - vep/$type.only_intron.txt >vep/$type.noncoding_notintron_affecting.txt
  awk '{print $1"\tfull_overlap"}' vep/$type.fulloverlap_affecting.txt |cat - <(awk '{print $1"\tpartial_CDS"}' vep/$type.partialCDS_affecting.txt) |cat - <(awk '{print $1"\tintronic"}' vep/$type.only_intron.txt) |cat - <(awk '{print $1"\tnoncoding"}' vep/$type.noncoding_notintron_affecting.txt) >vep/$type.comb.detail.txt
done
