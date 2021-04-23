
rule drip_bigwig_to_bedgraph:
    conda:
        '../envs/ucsc.yml'
    input:
        'rawdata/DRIP/{sample}_{strand}.bw'
    output:
        'output/DRIP_bedlike/{sample}_{strand}.bedgraph'
    shell:'''
    bigWigToBedGraph {input} {output}
    '''


rule bedgraph_to_bed:
    # input bedgraph with chr, start, end and score and return bed
    # with input, start, end, name (sample plus line number) score and strand
    input:
        'output/DRIP_bedlike/{sample}_{strand}.bedgraph'
    output:
        'output/DRIP_bedlike{sample}_{strand}.bed'
    shell:'''
    awk '{{print $1 "\t" $2 "\t" $3 "\t" NR "_{sample}" "\t" $4 "\t" $5}}'
    '''


rule sort_DRIP_bedlike:
    input:
        'output/DRIP_bedlike{sample}_{strand}.bed'
    output:
        'output/DRIP_bedlike/{sample}_{strand}.sorted.bed'
    shell:'''
    sort -k1,1n -k2,2n {input} > {output}
    '''


