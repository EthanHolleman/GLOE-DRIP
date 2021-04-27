

treatment_control_dict = {
    'siRNA LIG1': 'siRNA control',
    'LIG3-/-:mL3 + siRNA LIG1': 'LIG3-/-:mL3 + siRNA control'
}

def treatment_control_pairs(*args, **kwargs):
    '''Match up the treatment samples with their control samples so can be
    used in an expand statement later on which directs which samples get peak
    called together. Ouputs a list containing treatment control pairs, within
    each pair are two lists the first containing the treatment samples and
    the second containing the control samples. There are two samples per
    condition.
    '''
    pairs = []
    for treatment, control in treatment_control_dict.items():
        treatment_reps = list(GLOE_SAMPLES.loc[GLOE_SAMPLES['modification'] == treatment]['Sample Name'])
        control_reps = list(GLOE_SAMPLES.loc[GLOE_SAMPLES['modification'] == control]['Sample Name'])
        pairs.append([treatment_reps, control_reps])

    return pairs


GLOE_TRT_CNTRL_PAIRS = treatment_control_pairs()
STRANDS = ['fwd', 'rev']



rule seperate_forward_gloe_strands:
    input:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.sorted.bed'
    params:
        out_dir='output/call_gloe_peaks/strand'
    output:
        fwd='output/call_gloe_peaks/strand/{rep_a}_{rep_b}.fwd.sorted.bed'
    shell:'''
    mkdir -p {params.out_dir}
    grep "+" {input} > {output.fwd}
    '''


rule seperate_reverse_gloe_strands:
    input:
        'output/merged_gloe_replicates/{rep_a}_{rep_b}.sorted.bed'
    params:
        out_dir='output/call_gloe_peaks/strand'
    output:
        rev='output/call_gloe_peaks/strand/{rep_a}_{rep_b}.rev.sorted.bed'
    shell:'''
    mdkir {params.out_dir}
    grep "-" {input} > {output.rev}
    '''


rule call_peaks_treatment_as_treatment:
    conda:
        '../envs/macs2.yml'
    input:
        treatment='output/call_gloe_peaks/strand/{treat_rep_a}_{treat_rep_b}.{strand}.sorted.bed',
        control='output/call_gloe_peaks/strand/{control_rep_a}_{control_rep_b}.{strand}.sorted.bed'
    output:
        'output/call_gloe_peaks/TAT/treatment={treat_rep_a}_{treat_rep_b}.control={control_rep_a}_{control_rep_b}.{strand}._summits.bed'
    params:
        experiment_name='treatment={treat_rep_a}_{treat_rep_b}.control={control_rep_a}_{control_rep_b}.{strand}',
          out_dir='output/call_gloe_peaks/TAT'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.treatment} -c {input.control} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''


rule call_peaks_control_vs_treatment:
    conda:
        '../envs/macs2.yml'
    input:
        treatment='output/call_gloe_peaks/strand/{treat_rep_a}_{treat_rep_b}.{strand}.sorted.bed',
        control='output/call_gloe_peaks/strand/{control_rep_a}_{control_rep_b}.{strand}.sorted.bed'
    output:
        'output/call_gloe_peaks/CAT/treatment={control_rep_a}_{control_rep_b}.control={treat_rep_a}_{treat_rep_b}.{strand}._summits.bed'
    params:
        experiment_name='treatment={control_rep_a}_{control_rep_b}.control={treat_rep_a}_{treat_rep_b}.{strand}',
        out_dir='output/call_gloe_peaks/CAT'
    shell:'''
    mkdir -p {params.out_dir}
    macs2 callpeak -t {input.control} -c {input.treatment} -n {params.experiment_name} \
    --outdir {params.out_dir} -m 5 50 -g 1.20E+07 --format BED --extsize 1 \
    --nomodel --shift 0 --keep-dup all
    '''


rule call_all_peaks_treatment_as_treatment:
    input:
        expand(
            expand(
                'output/call_gloe_peaks/TAT/treatment={rep[0][0]}_{rep[0][1]}.control={rep[1][0]}_{rep[1][1]}.{strand}._summits.bed',
                 zip, rep=GLOE_TRT_CNTRL_PAIRS, allow_missing=True
            ),
           strand=STRANDS
        )
    output:
        'output/call_gloe_peaks/treatment_as_treatment.done'


rule call_all_peaks_control_as_treatment:
    input:
        expand(
            expand(
                'output/call_gloe_peaks/CAT/treatment={rep[1][0]}_{rep[1][1]}.control={rep[0][0]}_{rep[0][1]}.{strand}._summits.bed',
                zip, rep=GLOE_TRT_CNTRL_PAIRS, allow_missing=True
                ),
            strand=STRANDS
        )
    output:
        'output/call_gloe_peaks/control_as_treatment.done'
    shell:'''
    touch {output}
    '''


rule call_all_gloe_peaks:
    input:
        tat =  'output/call_gloe_peaks/treatment_as_treatment.done',
        cat = 'output/call_gloe_peaks/control_as_treatment.done'
    output:
        'output/call_gloe_peaks/all_gloe_peaks.done'
    shell:'''
    touch {output}
    '''
    



