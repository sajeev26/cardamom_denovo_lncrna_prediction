rule merge:
    input:
        cpc2="results/CPC2/cpc2.txt",
        plek="results/PLEK/plek"
    output:
        "results/lncRNA/lncRNA.csv"
    conda:
        "conda.yaml"
    shell:
        """
        mkdir -p results/lncRNA
        python script/merge.py \
            --cpc2 {input.cpc2} \
            --plek {input.plek} \
            --out {output}
        """
        
rule filter:
    input:
        "results/lncRNA/lncRNA.csv"
    output:
        "results/lncRNA/confident_lncRNA.csv"
    conda:
        "conda.yaml"
    shell:
        """
        mkdir -p results/lncRNA
        python script/filter.py \
            --in {input} \
            --out {output}
        """


rule extract:
    input:
        fasta="results/cdhit/transcripts.cdhit.fa",
        lnc="results/lncRNA/confident_lncRNA.csv"
    output:
        "results/lncRNA/confident_lncRNA.fasta"
    conda:
        "conda.yaml"
    shell:
        """
        python script/extract.py \
            --fasta {input.fasta} \
            --lnc {input.lnc} \
            --out {output}
        """

