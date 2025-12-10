# rules/qc.smk

rule fastqc_raw:
    input:
        r1="results/raw/{id}_1.fastq.gz",
        r2="results/raw/{id}_2.fastq.gz"
    output:
        html1="results/fastqc/raw/{id}_1_fastqc.html",
        html2="results/fastqc/raw/{id}_2_fastqc.html"
    conda:
    	"conda.yaml"
    shell:
        "fastqc {input.r1} {input.r2} -o results/fastqc/raw"

rule fastp:
    input:
        r1="results/raw/{id}_1.fastq.gz",
        r2="results/raw/{id}_2.fastq.gz"
    output:
        r1_trim="results/trimmed/{id}_1.trim.fastq.gz",
        r2_trim="results/trimmed/{id}_2.trim.fastq.gz",
        html="results/fastp/{id}_report.html",
        json="results/fastp/{id}_report.json"
    threads: 8
    conda:
    	"conda.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_trim} -O {output.r2_trim} \
              -w {threads} --detect_adapter_for_pe \
              --html {output.html} --json {output.json}
        """

rule fastqc_trimmed:
    input:
        r1="results/trimmed/{id}_1.trim.fastq.gz",
        r2="results/trimmed/{id}_2.trim.fastq.gz"
    output:
        html1="results/fastqc/trimmed/{id}_1_fastqc.html",
        html2="results/fastqc/trimmed/{id}_2_fastqc.html"
    conda:
    	"conda.yaml"
    shell:
        "fastqc {input.r1} {input.r2} -o results/fastqc/trimmed"
