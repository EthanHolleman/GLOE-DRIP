


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


rule metagene_promotors_plot:
    conda:
        '../envs/R.yml'
    input:
        bam_files=expand(
            'output/gloe_metaplots/{sample}.processed.direct.trimmed.sorted.bam',
            sample=SAMPLE_NAMES),
        indicies=expand(
            'output/gloe_metaplots/{sample}.processed.direct.trimmed.sorted.bai',
            sample=SAMPLE_NAMES)
    output:
        'output/gloe_metaplots/plots/promotor_metaplot.png'
    params:
        sample_names = SAMPLE_NAMES
    shell:'''
    Rscript scripts/metaplot_gloe_replicates.R {output} {input.bam_files} {input.sample_names}
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


