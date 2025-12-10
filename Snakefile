# Snakefile

include: "rules/download.smk"
include: "rules/qc.smk"
include: "rules/spades.smk"
include: "rules/quant.smk"
include: "rules/assembly_qc.smk"
include: "rules/cp.smk"
include: "rules/filter.smk"

configfile: "config/config.yaml"
SAMPLES = [line.strip() for line in open("data/sra.txt") if line.strip()]

rule all:
    input:
        "results/rnaspades/transcripts.fasta",
        "results/cdhit/transcripts.cdhit.fa",
        expand("results/salmon/{id}/quant.sf", id=SAMPLES),
        "results/busco/run_busco",
        "results/CPC2/cpc2.txt",
        "results/PLEK/plek",
        "results/lncRNA/lncRNA.csv",
        "results/lncRNA/confident_lncRNA.csv",
        "results/lncRNA/confident_lncRNA.fasta"
