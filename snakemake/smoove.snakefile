REF = "data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
GFF = "data/raw_data/reference/GCF_000002285.3_CanFam3.1_genomic.gff",
SAMPLES = ["AfGWo", "BlBJa"]   # specify samples within this pipeline

rule all:
    input:
        #expand("data/bam/{prefix}.md.bam.bai", prefix=SAMPLES)
        "smoove/genotyped/BlBJa.joint-smoove.genotyped.vcf.gz"
        "smoove/genotyped/AfGWo.joint-smoove.genotyped.vcf.gz"

#rule call:
#    input:
#        #bam="sorted_reads/{wolf}.bam"
#        cram="data/raw_data/cram/{wolf}.cram"
#    output:
#        "out/smoove/results/called/{wolf}-smoove.genotyped.vcf.gz"
#    log:
#        "logs/call/{wolf}.log"
#    singularity:
#        "smoove_latest.sif"
#    shell:
        #"smoove call --outdir smoove/results/ --exclude {EXCLUDE} –name {wildcards.wolf} --fasta {REF} -p 1 --genotype {input.bam}"
#        "smoove call --outdir out/smoove/results/called/ --exclude {EXCLUDE} --name {wildcards.wolf} --fasta {REF} -p 8 --genotype {input.cram}"
#rule merge:
#    input:
#        vcf=expand("out/smoove/results/called/{wolf}-smoove.genotyped.vcf.gz",
#        wolf=WOLVES)
#    output:
#        "out/smoove/results/called/merged.sites.vcf.gz"
#    log:
#        "logs/merge/merged.sites.log"
#    singularity:
#        "smoove_latest.sif"
#    shell:
#        "smoove merge --name merged -f {REF} --outdir out/smoove/results/called/ {input.vcf}"

rule genotype:
    input:
        merge="../../lars/thesis/data/smoove/called/merged.sites.vcf.gz",
        bam="data/bam/{prefix}.md.bam"
        #cram="data/raw_data/cram/{wolf}.cram"
    output:
        "smoove/genotyped/{prefix}.joint-smoove.genotyped.vcf.gz"
    log:
        "logs/genotype/{prefix}.log"
    singularity:
        "smoove_latest.sif"
    shell:
            "smoove genotype -d -p 8 --name {wildcards.prefix}.joint --outdir smoove/genotyped/ --fasta {REF} --vcf {input.merge} {input.bam}"

#rule paste:
#    input:
#        vcf=expand("out/smoove/results/genotyped/{wolf}.joint-smoove.genotyped.vcf.gz",
#        wolf=WOLVES)
#    output:
#        "cohort.smoove.square.vcf.gz"
#    log:
#        "logs/paste/cohort.log"
#    singularity:
#        "smoove_latest.sif"
#    shell:
#        "smoove paste --name cohort {input.vcf}"
#
#rule annotate:
#    input:
#        "cohort.smoove.square.vcf.gz"
#    output:
#        "out/smoove/results/annotated/cohort.smoove.square.anno.vcf.gz"
#    log:
#        "logs/annotate/cohort.smoove.anno.log"
#    singularity:
#        "smoove_latest.sif"
#    shell:
#        "smoove annotate --gff {GFF} {input} | bgzip -c > {output}"
#
