import pandas as pd

GLOE_SAMPLES = pd.read_csv(
    'samples/GLOE_samples.csv', sep=','
).set_index('Sample Name', drop=False)
SAMPLE_NAMES=GLOE_SAMPLES['Sample Name']

include: 'rules/download.smk'
include: 'rules/trim_gloe_seq.smk'
include: 'rules/map_gloe_reads.smk'
include: 'rules/process_mapped_gloe_reads.smk'

wildcard_constraints:
    sample='\w+'


MAP_GLOE_READS = 'output/map_glow_reads/map_all_reads.done'
PROCESS_GLOE_READS = 'output/process_alignment/process_all_alignments.done'
MERGE_GLOE_REPLICATES = 'output/merged_gloe_replicates/merge_all_gloe_replicates.done'

rule all:
    input:
        MAP_GLOE_READS,
        PROCESS_GLOE_READS, 
        MERGE_GLOE_REPLICATES
       
        