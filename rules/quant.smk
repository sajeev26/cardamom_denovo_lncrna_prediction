#rules/quant.smk

rule salmon_index:
    input:
        "results/cdhit/transcripts.cdhit.fa"
    output:
        directory("results/salmon/index")
    conda:
    	"conda.yaml"
    shell:
        "salmon index -t {input} -i {output} -k 31"

rule salmon_quant:
    input:
        index="results/salmon/index",
        r1="results/trimmed/{id}_1.trim.fastq.gz",
        r2="results/trimmed/{id}_2.trim.fastq.gz"
    output:
        "results/salmon/{id}/quant.sf"
    conda:
    	"conda.yaml"
    shell:
        """
        salmon quant -i {input.index} \
                     -l A -1 {input.r1} -2 {input.r2} \
                     -p 8 -o results/salmon/{wildcards.id}
        """
