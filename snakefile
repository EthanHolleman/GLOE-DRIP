import pandas as pd

GLOE_SAMPLES = pd.read_csv(
    'samples/GLOE_samples.csv', sep=','
).set_index('Sample Name', drop=False)
SAMPLE_NAMES=GLOE_SAMPLES['Sample Name']

DRIP_SAMPLES = pd.read_csv(
    'samples/DRIP_samples.tsv', sep='\t'
).set_index('sample', drop=False)




include: 'rules/download.smk'
include: 'rules/trim_gloe_seq.smk'
include: 'rules/map_gloe_reads.smk'
include: 'rules/process_mapped_gloe_reads.smk'
include: 'rules/merge_gloe_replicates.smk'
include: 'rules/call_gloe_peaks.smk'
include: 'rules/drip_to_bedlike.smk'
include: 'rules/metaplot_gloe_replicates.smk'

wildcard_constraints:
    sample='\w+',
    rep_a='\w+',
    rep_b='\w+',
    strand='\w+',
    treat_rep_a='\w+',
    treat_rep_b='\w+',
    control_rep_a='\w+',
    control_rep_b='\w+'

DOWNLOAD_DRIP_DATA = 'rawdata/DRIP/download_all_drip.done'
DRIP_TO_BED = 'output/DRIP_bedlike/all_drip_to_bedlike.done'
MAP_GLOE_READS = 'output/map_glow_reads/map_all_reads.done'
PROCESS_GLOE_READS = 'output/process_alignment/process_all_alignments.done'
MERGE_GLOE_REPLICATES = 'output/merged_gloe_replicates/merge_all_gloe_replicates.done'
CALL_GLOE_PEAKS = 'output/call_gloe_peaks/all_gloe_peaks.done'

rule all:
    input:
        DOWNLOAD_DRIP_DATA,
        DRIP_TO_BED,
        #MAP_GLOE_READS,
        #PROCESS_GLOE_READS, 
        #MERGE_GLOE_REPLICATES,
        #CALL_GLOE_PEAKS,
        'output/gloe_metaplots/index_all_bams.done'
       
        