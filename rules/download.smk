# rules/download.smk


rule download_sra:
    input:
        "data/sra.txt"
    output:
        temp("results/raw/{id}_1.fastq.gz"),
        temp("results/raw/{id}_2.fastq.gz")
    threads: 4
    conda:
    	"conda.yaml"
    shell:
        """
        fasterq-dump --split-files -O results/raw {wildcards.id}
        gzip results/raw/{wildcards.id}_1.fastq
        gzip results/raw/{wildcards.id}_2.fastq
        """
