

rule truncate_drip_bedg_fwd:
    input:
        'output/call_drip_peaks/{sample}_fwd/{sample}_fwd_summits.bed'
    output:
        'output/truncated_drip/{sample}_fwd.truncated.bed'
    params:
        extra_bases=10
    shell:'''
    awk '{{OFS="\t" print $1,$2,$2+{params.extra_bases},$4,$5,$6}}'
    '''


rule truncate_drip_bed_rev:
    # truncate from end since initiation site is really at the end of
    # the sequence since flipped
    input:
        'output/call_drip_peaks/{sample}_rev/{sample}_rev_summits.bed'
    output:
        'output/truncated_drip/{sample}_rev.truncated.bed'
    params:
        extra_bases=10
    shell:'''
    awk '{{OFS="\t" print $1,$3-{params.extra_bases},$3,$4,$5,$6}}'
    '''


rule combine_truncated_strands:
    input:
        rev='output/truncated_drip/{sample}_rev.truncated.bed',
        fwd='output/truncated_drip/{sample}_fwd.truncated.bed'
    output:
        'output/truncated_drip/{sample}_all.truncated.bed'
    shell:'''
    cat {input.fwd} {input.rev} > {output}
    '''


rule sort_combined_truncated_strands:
    input:
        'output/truncated_drip/{sample}_all.truncated.bed'
    output:
        'output/truncated_drip/{sample}_all.truncated.sorted.bed'
    shell:'''
    sort -k1,1n -k2,2n {input} > {output}
    '''



