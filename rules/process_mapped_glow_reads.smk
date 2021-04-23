
rule sam_to_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/bowtie2/{sample}.sam'
    output:
        'output/{sample}/process_alignment/{sample}.bam'
    shell:'''
    samtools view -bhSu -o {output} {input}
    '''

rule sort_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/process_alignment/{sample}.bam'
    output:
        temp('output/{sample}/process_alignment/{sample}.sorted.bam')
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
        'output/{sample}/process_alignment/{sample}.sorted.bam'
    output:
        'output/{sample}/process_alignment/trim_low_qual/all.bam'
    shell:'''
    samtools view -q 30 -bhu -o {output} {input}
    '''


rule sort_trimmed_bam:
    conda: 
        '../envs/samtools.yml'
    input:
        'output/{sample}/process_alignment/trim_low_qual/{region}.bam'
    output:
        'output/{sample}/process_alignment/sorted/trim.{region}.bam'
    threads: 16
    params:
        sort='output/{sample}/alignment/{sample}.{region}.trim.temp'
    shell:'''
    mkdir --parents {params.sort}
    samtools sort -O bam -T {params.sort} --threads {threads} -o {output} {input}
    rm -r {params.sort}
    '''


rule bam_to_bed:
    conda: 
        '../envs/bedtools.yml'
    input:
        'output/{sample}/process_alignment/sorted/trim.{region}.bam'
    output:
        temp('output/{sample}/process_alignment/bed/sorted.trim.{region}.bed')
    shell:'''
    bedtools bamtobed -i {input} > {output}
    '''


rule direct_mode:
    input:
        'output/{sample}/process_alignment/bed/sorted.trim.{region}.bed'
    output:
        sites=temp('output/{sample}/reorient_alignments/direct/trim.{region}.bed')
    params:
        index_dir='output/{sample}/direct'
    shell:'''
    module load perl
    mkdir -p {params.index_dir}
    perl scripts/direct_mode.pl {input} > {output.sites}
    '''


rule get_second_column:
    input:
        'output/{sample}/reorient_alignments/{mode}/trim.{region}.bed'
    output:
        temp('output/{sample}/reorient_alignments/{mode}/sorted.col2.trim.{region}.bed')
    shell:"""
    awk '($2 >= 0)' {input} > {output}
    """


rule perl_mode_big_awk:
    input:
        'output/{sample}/reorient_alignments/{mode}/sorted.col2.trim.{region}.bed'
    output:
        'output/{sample}/reorient_alignments/{mode}/bigawk.sorted.trim.{region}.bed'
    shell:"""
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}}'  {input} > {output}
    """

