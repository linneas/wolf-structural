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

args = parser.parse_args()

f = open(args.contigs, 'r')
contigs = json.load(f)
f.close()

vcf = VariantFile(args.vcf)
samples = list(vcf.header.samples)

for sv in vcf:

    if sv.contig not in contigs:
        continue

    if sv.info["SVTYPE"] != 'DEL':
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

        file_name = sv.info['SVTYPE'] + '_' \
                    + sv.contig + '_' \
                    + str(sv.start) + '_' \
                    + str(sv.stop) + '_' \
                    + '_'.join([str(x) for x in sample_idxs])

        print(file_name, sv.contig, str(sv.start), str(sv.stop), sv.info['SVTYPE'], sep='\t')
