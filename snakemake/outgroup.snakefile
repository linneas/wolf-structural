# Written by LS
SAMPLES = ["AfGWo", "BlBJa"]   # specify samples within this pipeline

# The following rules will not be executed as jobs
localrules: all, makeGVCFList,

# Without specifying a file, Snakemake runs the first rule by default.
# Therefore, specify a rule that has the final output as input -> this will
# run the entire pipeline!
rule all:
    input:
        expand("flagstats/{prefix}.stats", prefix=SAMPLES),
        expand("data/bam/{prefix}.md.bam.bai", prefix=SAMPLES)

# Run a rule from a separate file
# include: "mapping.snakefile"

# the r in the rg parameter makes avoids replacing "\t" for tabs in the command
rule bwa:
    input:
        ref="data/raw_data/reference/Canis_familiaris.CanFam3.1.dna.toplevel.fa",
        r1="data/raw_data/fastq/{prefix}_1.fastq.gz",
        r2="data/raw_data/fastq/{prefix}_2.fastq.gz"
    output:
        "data/bam/{prefix}.bam"
    log:
        "logs/bwa/{prefix}.log"
    params:
        rg=r"@RG\tID:{prefix}\tSM:{prefix}\tLB:{prefix}\tPL:Illumina"
    #    tmp=r"$SNIC_TMP/{prefix}"
    threads: 20
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} -M {input.ref} {input.r1} {input.r2} \
         |samtools sort -m 6G -@{threads} -T $SNIC_TMP/{wildcards.prefix} - >{output}.tmp) 2>{log} \
         && mv {output}.tmp {output}"

rule index_bam:
    input:
        "data/bam/{prefix}.bam"
    output:
        "data/bam/{prefix}.bam.bai"
    log:
        "logs/samtools/index.{prefix}.log"
    threads: 1
    shell:
        "samtools index {input} {output}"

rule flagstat:
    input:
        "data/bam/{prefix}.bam"
    output:
        "flagstats/{prefix}.stats"
    log:
        "logs/samtools/flagstats.{prefix}.log"
    threads: 1
    shell:
        "(samtools flagstat {input} >{output}) 2>{log}"

rule mkdupl:
    input:
        "data/bam/{prefix}.bam"
    output:
        bam="data/bam/{prefix}.md.bam",
        metrics="data/bam/{prefix}.metrics"
    log:
        "logs/picard/{prefix}.mkdupl.log"
    params:
        picard="/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar",
        mem="60g"
    threads: 10
    shell:
        "(java -Xmx{params.mem} -jar {params.picard} MarkDuplicates INPUT={input} \
         METRICS_FILE={output.metrics} TMP_DIR=$SNIC_TMP ASSUME_SORTED=true \
         VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT={output.bam}.tmp) 2>{log} \
         && mv {output.bam}.tmp {output.bam} && mv {output.bam}.tmp.bai {output.bam}.bai"
