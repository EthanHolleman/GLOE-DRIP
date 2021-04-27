


rule bed_to_bam:
    # convert processed bed files back to bam
    conda:
        '../envs/bedtools.yml'
    input:
        bed='output/process_alignment/{sample}.processed.direct.bed',
        genome='rawdata/hg19/hg19.chrom.sizes'
    output:
        'output/gloe_metaplots/{sample}.processed.direct.bam'
    params:
        out_dir='output/gloe_metaplots'
    shell:'''
    mkdir -p {params.out_dir}
    bedtools bedtobam -i {input.bed} -g {input.genome} > {output}
    '''


rule sort_bam_for_index:
    conda:
        '../envs/samtools.yml'
    input:
        'output/gloe_metaplots/{sample}.processed.direct.bam'
    output:
        'output/gloe_metaplots/{sample}.processed.direct.sorted.bam'
    shell:'''
    samtools sort {input} -o {output}
    '''


rule index_bam:
    conda:
        '../envs/samtools.yml'
    input:
        'output/gloe_metaplots/{sample}.processed.direct.sorted.bam'
    output:
        'output/gloe_metaplots/{sample}.processed.direct.sorted.bai'
    shell:'''
    samtools index {input} {output}
    '''

rule make_gloe_design_file:
    input:
        'samples/GLOE_samples.csv'
    output:
        'output/gloe_metaplots/design/gloe_design.tsv'
    params:
        out_dir = 'output/gloe_metaplots/design'
    shell:'''
    mkdir -p {params.out_dir}
    python scripts/design.py {input} {output}
    '''


rule metagene_promotors_plot:
    conda:
        '../envs/R.yml'
    input:
        bam_files=expand(
            'output/gloe_metaplots/{sample}.processed.direct.sorted.bam',
            sample=SAMPLE_NAMES),
        indicies=expand(
            'output/gloe_metaplots/{sample}.processed.direct.sorted.bai',
            sample=SAMPLE_NAMES),
        design = 'output/gloe_metaplots/design/gloe_design.tsv'
    output:
        'output/gloe_metaplots/plots/promotor_metaplot.png'
    params:
        sample_names = lambda wildcards: list(SAMPLE_NAMES)
    shell:'''
    Rscript scripts/metagene.R {output} {input.design} {input.bam_files} {params.sample_names}
    '''
        
        
rule index_all_bams:
    input:
        expand('output/gloe_metaplots/{sample}.processed.direct.sorted.bai',
        sample=SAMPLE_NAMES)
    output:
        'output/gloe_metaplots/index_all_bams.done'
    shell:'''
    touch {output}
    '''


