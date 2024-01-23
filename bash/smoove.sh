#!/bin/bash

# Code by Lars Huson, adapted from Bertolotti et al 2020

# Software versions
# snakemake/7.18.2

# create flowchart of snakemake pipeline
snakemake --snakefile snakemake/Snakefile --dag | dot -Tsvg > snake_dag.svg

# prep slurm dir
mkdir logs/slurm
mkdir smoove/results/annotated
# run SV caller in snakemake
snakemake -s snakemake/Snakefile --use-singularity --nolock -j 64 \
        --cluster "sbatch -A p2018002 -p {cluster.partition} -n {cluster.n} -t {cluster.time} -J {cluster.name} --parsable -e {cluster.error}.slurm -o {cluster.output}.slurm --mail-user {cluster.email} --mail-type {cluster.mailtype}" \
        --cluster-config snakemake/config_snake_cluster.json
#NOTE: Using the current config file, jobs for the 'call' step in snakemake will be run on 20 cores to accommodate for wolf 56-D-07-16.
#       All other jobs can be run successfully on 12 cores or less.
