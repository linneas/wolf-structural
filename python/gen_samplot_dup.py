# Script provided by Gabriel David

import argparse
from pysam import VariantFile
import json
from collections import namedtuple
import random

Genotypes = namedtuple('Genotypes', 'REF HET ALT UNK')

parser = argparse.ArgumentParser()

parser.add_argument('--vcf',
                    dest='vcf',
                    required=True)

parser.add_argument('--contigs',
                    dest='contigs',
                    required=True)

parser.add_argument('--bams',
                    dest='bams',
                    required=True)

parser.add_argument('--outdir',
                    dest='outdir',
                    required=True)

args = parser.parse_args()

f = open(args.contigs, 'r')
contigs = json.load(f)
f.close()

f = open(args.bams, 'r')
bams = json.load(f)
f.close()

vcf = VariantFile(args.vcf)
samples = list(vcf.header.samples)

for sv in vcf:

    if sv.contig not in contigs:
        continue

    if sv.info["SVTYPE"] != 'DUP':
        continue

    gts = Genotypes(REF=[], HET=[], ALT=[], UNK=[])
    for s in samples:
        if sv.samples[s]['GT'] == (0,0):
            gts.REF.append(s)
        elif sv.samples[s]['GT'] == (0,1):
            gts.HET.append(s)
        elif sv.samples[s]['GT'] == (1,1):
            gts.ALT.append(s)
        elif sv.samples[s]['GT'] == (None,None):
            gts.UNK.append(s)



    if len(gts.REF) >= 2 and len(gts.HET) >= 2 and len(gts.ALT) >= 2:
        rand_refs = random.sample(gts.REF, k=2)
        rand_hets = random.sample(gts.HET, k=2)
        rand_alts = random.sample(gts.ALT, k=2)

        sample_idxs = []
        for s in rand_refs + rand_hets + rand_alts:
            sample_idxs.append(samples.index(s))

        bam_paths = []
        for s in rand_refs + rand_hets + rand_alts:
            bam_paths.append(bams[s])

        file_name = sv.info['SVTYPE'] + '_' \
                    + sv.contig + '_' \
                    + str(sv.start) + '_' \
                    + str(sv.stop) + '_' \
                    + '_'.join([str(x) for x in sample_idxs]) \
                    + '.png'

        print('samplot plot -r data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa' \
              + ' -n ' + ' '.join(rand_refs + rand_hets + rand_alts) \
              + ' -b ' + ' '.join(bam_paths) \
              + ' -o ' + 'imgs_DUPs/'+ file_name \
              + ' -c ' + sv.contig \
              + ' -s ' + str(sv.start) \
              + ' -e ' + str(sv.stop) \
              + ' -a ' \
              + ' -t ' + sv.info['SVTYPE'])
