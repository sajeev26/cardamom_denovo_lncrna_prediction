# rules/spades.smk

SAMPLES = [line.strip() for line in open("data/sra.txt") if line.strip()]


rule rnaspades_assembly:
    input:
        r1 = expand("results/trimmed/{id}_1.trim.fastq.gz", id=SAMPLES),
        r2 = expand("results/trimmed/{id}_2.trim.fastq.gz", id=SAMPLES)
    output:
        "results/rnaspades/transcripts.fasta"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        mkdir -p results/rnaspades

        # CONCATENATE ALL R1 AND R2 READS
        cat {input.r1} > results/rnaspades/all_R1.fastq.gz
        cat {input.r2} > results/rnaspades/all_R2.fastq.gz

        rnaspades.py \
            -1 results/rnaspades/all_R1.fastq.gz \
            -2 results/rnaspades/all_R2.fastq.gz \
            -o results/rnaspades/rnaspades_out \
            --threads {threads}

        cp results/rnaspades/rnaspades_out/transcripts.fasta {output}
        """


rule cdhit_reduce:
    input:
        "results/rnaspades/transcripts.fasta"
    output:
        "results/cdhit/transcripts.cdhit.fa"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        mkdir -p results/cdhit
        cd-hit-est \
            -i {input} \
            -o {output} \
            -c 0.95 \
            -T {threads} \
            -M 0
        """


