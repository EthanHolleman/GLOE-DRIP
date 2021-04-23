

rule download_all_gloe_samples:
    conda:
        '../envs/sra-toolkit.yml'
    output:
       temp('rawdata/GLOE-seq/{sample_name}.sra')
    params:
        sample_name = lambda wildcards: wildcards['sample_name'],
        sra_accession = lambda wildcards: GLOE_SAMPLES.loc[wildcards.sample_name]['Run']
    shell:'''
    rm -rf rawdata/GLOE-seq/{params.sample_name}.sra.lock
    prefetch {params.sra_accession} --output-file {output}
    '''


rule dump_gloe_fastq:
    input:
        'rawdata/GLOE-seq/{sample}.sra'
    output:
        'rawdata/GLOE-seq/{sample}.fastq.gz'
    shell:'''
    fastq-dump -Z {input} | gzip > {output}
    '''


# Download primers
rule download_primer_file:
    output:
        'rawdata/primers/TruSeq3-SE.fa'
    shell:'''
    curl https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa \
    -o {output}
    '''

rule download_hg19_chr_sizes:
    output:
        'rawdata/hg19/hg19.chrom.sizes'
    shell:'''
    curl -L http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes -o {output}
    '''