
rule drip_bigwig_to_bedgraph:
    conda:
        '../envs/ucsc.yml'
    input:
        'rawdata/DRIP/{sample_strand}.bw'
    output:
        'output/DRIP_bedlike/{sample_strand}.bedgraph'
    shell:'''
    bigWigToBedGraph {input} {output}
    '''


rule bedgraph_to_bed:
    # input bedgraph with chr, start, end and score and return bed
    # with input, start, end, name (sample plus line number) score and strand
    input:
        'output/DRIP_bedlike/{sample_name}_{strand}.bedgraph'
    output:
        'output/DRIP_bedlike{sample_name}_{strand}.bed'
    params:
        strand = lambda wildcards: '+' if wildcards['strand'] == 'fwd' else '-',
        sample_name = lambda wildcards: wildcards['sample_name']
    shell:'''
    awk '{{print $1 "\t" $2 "\t" $3 "\t" NR "_{params.sample_name}" "\t" $4 "\t" "{params.strand}"}}' {input} > {output}
    '''


rule sort_DRIP_bedlike:
    input:
        'output/DRIP_bedlike{sample}_{strand}.bed'
    output:
        'output/DRIP_bedlike/{sample}_{strand}.sorted.bed'
    shell:'''
    sort -k1,1n -k2,2n {input} > {output}
    '''


rule all_drip_to_bedlike:
    input:
        expand('output/DRIP_bedlike/{sample}.sorted.bed', 
        sample=DRIP_SAMPLES['sample']
        )
    output:
        'output/DRIP_bedlike/all_drip_to_bedlike.done'
    shell:'''
    touch {output}
    '''


