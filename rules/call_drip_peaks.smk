

rule call_drip_peaks:
    conda:
        '../envs/macs2.yml'
    input:
        'output/DRIP_bedgraph/{sample}_{strand}.sorted.bedgraph'
    output:
        'output/call_drip_peaks/{sample}_{strand}/{sample}_{strand}_summits.bed
    params:
        experiment_name='{sample}_{strand}'
        out_dir='output/call_drip_peaks/{sample}_{strand}'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 peakcall -t {input} -n {params.experiment_name} --outdir {params.outdir} \
    -m 5 -g 1.20E+07 --format BED --nomodel --shift 0
    '''

rule call_all_drip_peaks:
    input:
        expand(
            'output/call_drip_peaks/{sample}_{strand}/{sample}_{strand}_summits.bed',
            sample
        )


    