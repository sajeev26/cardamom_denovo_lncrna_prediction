rule busco:
    input:
        "results/cdhit/transcripts.cdhit.fa"
    output:
        directory("results/busco/run_busco")
    conda:
        "conda.yaml"
    shell:
        """
        mkdir -p results/busco

        busco \
            -i {input} \
            -l embryophyta_odb10 \
            -m transcriptome \
            -o run_busco \
            -c 8
        """

