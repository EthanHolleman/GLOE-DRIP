

rule concatenate_gloe_replicates:
    input:
        rep_a='',
        rep_b=''
    output:     
        'output/merged_gloe_replicates/{sample_a}_{sample_b}.bed'
    shell:'''
    cat {input.rep_a} {input.rep_b} > {output}
    '''


rule merge_gloe_replicates:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/merged_gloe_replicates/{sample_a}_{sample_b}.bed'
    output:
        'output/merged_gloe_replicates/{sample_a}_{sample_b}.bed'
    shell:'''
    bedtools merge 
