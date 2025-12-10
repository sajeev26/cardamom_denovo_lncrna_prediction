# workflow/rules/coding_potential.smk

rule cpc2:
    input:
        "results/cdhit/transcripts.cdhit.fa"
    output:
        "results/CPC2/cpc2.txt"
    conda:
    	"conda.yaml"
    shell:
        """
        python CPC2/bin/CPC2.py -i {input} -o results/CPC2/cpc2
        """

rule plek:
    input:
        "results/cdhit/transcripts.cdhit.fa"
    output:
        "results/PLEK/plek"
    conda:
    	"conda.yaml"
    shell:
        """
        python PLEK.1.2/PLEK.py -fasta {input} -out results/PLEK/plek -thread 8
        """



