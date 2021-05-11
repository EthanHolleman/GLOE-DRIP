# Read an SRA table with sample information and create a "design"
# table that can be used by the R metagene package for plotting
# by experimental group. In the sra file the experimental group
# should be under a column called "modification" and the sample
# name should be under "Sample Name" (case sensitive).

import pandas as pd 
import argparse
from pathlib import Path

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('sra', help='Path to sra table to read.')
    parser.add_argument('out', help='Path to output table.')
    parser.add_argument('--delim', default=',', help='Delimiter of sra file.')

    return parser.parse_args()


def read_sra(filepath, delim=',', **kwargs):
    return pd.read_csv(str(filepath), sep=delim, **kwargs)


def design_frame(sra_frame):
    groups = list(set(list(sra_frame['modification'])))
    groups = [g.replace(' ', '_') for g in groups]
    print(groups)
    group_matrix = {}
    for index, row in sra_frame.iterrows():
        sample_row = {g: 0 for g in groups}
        sample_row[str(row['modification']).replace(' ','_')] = 1
        group_matrix[str(row['Sample Name'])] = sample_row
    return pd.DataFrame(group_matrix).transpose()


def main():
    args = get_args()
    sra_df = read_sra(args.sra)
    design_df = design_frame(sra_df)
    design_df.to_csv(args.out, sep='\t')

if __name__ == '__main__':
    main()
    
        


    