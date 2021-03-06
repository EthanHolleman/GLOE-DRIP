
rule sam_to_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/map_glow_reads/{sample}.sam'
    output:
        temp('output/process_alignment/{sample}.bam')
    shell:'''
    mkdir -p output/process_alignment
    samtools view -bhSu -o {output} {input}
    '''

rule sort_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/process_alignment/{sample}.bam'
    output:
        temp('output/process_alignment/{sample}.sorted.bam')
    params:
        sort='output/{sample}_temp'
    threads: 16
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule trim_bam_low_qual_alignments_all:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/process_alignment/{sample}.sorted.bam'
    output:
        'output/process_alignment/{sample}.trimmed.bam'
    shell:'''
    samtools view -q 30 -bhu -o {output} {input}
    '''


rule sort_trimmed_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/process_alignment/{sample}.trimmed.bam'
    output:
        'output/process_alignment/{sample}.trimmed.sorted.bam'
    threads: 8
    params:
        sorting_temp='output/process_alignment/{sample}'
    shell:'''
    mkdir --parents {params.sorting_temp}
    samtools sort -O bam -T {params.sorting_temp} --threads {threads} -o {output} {input}
    rm -r {params.sorting_temp}
    '''


rule bam_to_bed:
    conda: 
        '../envs/bedtools.yml'
    input:
        'output/process_alignment/{sample}.trimmed.sorted.bam'
    output:
        temp('output/process_alignment/{sample}.bed')
    shell:'''
    bedtools bamtobed -i {input} > {output}
    '''


rule direct_mode:
    input:
        'output/process_alignment/{sample}.bed'
    output:
        sites='output/process_alignment/{sample}.direct.bed'
    shell:'''
    module load perl
    perl scripts/direct_mode.pl {input} > {output.sites}
    '''


rule get_second_column:
    input:
        'output/process_alignment/{sample}.direct.bed'
    output:
        temp('output/process_alignment/{sample}.nonzero.direct.bed')
    shell:"""
    awk '($2 >= 0)' {input} > {output}
    """


rule perl_mode_big_awk:
    input:
        'output/process_alignment/{sample}.nonzero.direct.bed'
    output:
        'output/process_alignment/{sample}.processed.direct.bed'
    shell:"""
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}}'  {input} > {output}
    """


rule sort_processed_bed:
    input:
        'output/process_alignment/{sample}.processed.direct.bed'
    output:
        'output/process_alignment/{sample}.processed.direct.sorted.bed'
    shell:'''
    sort -k 1,1 -k2,2n {input} > {output}
    '''


rule process_all_alignments:
    # ask for this file in the snakefile to run all rules in this file
    input:
        expand('output/process_alignment/{sample}.processed.direct.bed', 
        sample=GLOE_SAMPLES['Sample Name']
        )
    output:
        'output/process_alignment/process_all_alignments.done'
    shell:'''
    touch {output}
    '''
    

