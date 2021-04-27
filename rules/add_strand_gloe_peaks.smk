# Macs2 thinkss that there is no strandedness for peak calling
# so need to add it back in before concatenting

rule add_strand_to_summit_file:
    input:
        'output/call_gloe_peaks/CAT/treatment={sample_a}_{sample_b}.control={sample_c}_{sample_d}.{strand}._summits.bed'
    output:
        'output/stranded_gloe_summits/treatment={sample_a}_{sample_b}.control={sample_c}_{sample_d}.with.{strand}._summits.bed'
    params:
        out_dir = 'output/stranded_gloe_summits'
        strand=lambda wildcards: '+' if wildcards['strand'] == 'fwd' else '-'
    shell:'''
    mkdir -p {params.out_dir}
    awk '{{$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" {params.strand}}}' {inout} \
    > {output}
    '''

