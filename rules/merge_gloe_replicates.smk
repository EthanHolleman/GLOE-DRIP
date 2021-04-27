
def group_replicates(samples):
    all_samples = []
    for _, sample in samples.iterrows():
        modification = sample['modification']
        replicates = samples.loc[samples['modification'] == modification]
        replicates = tuple([rep['Sample Name'] for _, rep in replicates.iterrows()])
        all_samples.append(replicates)
    return all_samples


REPLICATES = group_replicates(GLOE_SAMPLES)


rule concatenate_gloe_replicates:
    input:
        rep_a='output/process_alignment/{rep_a}.processed.direct.trimmed.sorted.bed',
        rep_b='output/process_alignment/{rep_b}.processed.direct.trimmed.sorted.bed'
    output:     
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.bed'
    shell:'''
    mkdir -p output/merged_gloe_replicates
    cat {input.rep_a} {input.rep_b} > {output}
    '''

rule sort_concat_gloe_reps:
    input:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.bed'
    output:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.sorted.bed'
    shell:'''
   sort -k 1,1 -k2,2n {input} > {output}
    '''


rule merge_gloe_replicates:
    # strand specific merger 
    conda:
        '../envs/bedtools.yml'
    input:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.sorted.bed'
    output:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.merged.bed'
    shell:'''
    mkdir -p output/merged_gloe_replicates
    bedtools merge -s -c 5 -o average -i {input} > {output}
    '''


rule sort_merged_gloe_replicates:
    input:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.merged.bed'
    output:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.merged.sorted.bed'
    shell:'''
    sort -k1,1 -k2,2n {input} > {output}
    '''


rule merge_all_gloe_replicates:
    input:
        expand(
            'output/merged_gloe_replicates/{rep[0]}_{rep[1]}.merged.sorted.bed',
            zip, rep=REPLICATES
            )
    output:
        'output/merged_gloe_replicates/merge_all_gloe_replicates.done'
    
