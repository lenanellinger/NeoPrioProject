import os
import pandas as pd
import numpy as np
import sys
from optparse import OptionParser

def main(argv):
    usage = "usage: python neoantigen_prioritization.py "
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--pvacseq-filtered-output", action="store", dest="pvacseq_filtered_output", help="Filtered TSV file output of pVACseq")
    parser.add_option("--neofox-output", action="store", dest="neofox_output", help="TSV file output of neofox")
    (options, args) = parser.parse_args()
    
    pvacseq_df = pd.read_csv(options.pvacseq_filtered_output, sep="\t", header=0)
    neofox_df = pd.read_csv(options.neofox_output, sep="\t", header=0)
    
    if neofox_df.shape[0] != pvacseq_df.shape[0]:
        raise ValueError("Neofox output and pVACseq output should have the same amount of rows.")
    annotation_df = neofox_df.merge(right=pvacseq_df.loc[:, ['Chromosome', 'Start', 'Stop', 'Transcript', 'Variant Type', 'MHCflurry WT IC50 Score', 'MHCflurry MT IC50 Score', 'MHCflurry WT Percentile', 'MHCflurry MT Percentile', 
    'NetMHC WT IC50 Score', 'NetMHC MT IC50 Score', 'NetMHC WT Percentile', 'NetMHC MT Percentile', 
    'NetMHCpan WT IC50 Score', 'NetMHCpan MT IC50 Score', 'NetMHCpan WT Percentile',  'NetMHCpan MT Percentile', 
    'Predicted Stability', 'Half Life', 'Stability Rank', 'NetMHCstab allele']], on=['Chromosome', 'Start', 'Stop', 'Transcript'])
    
    # Prefilter (majority vote on binding affinities below 500)
    voting = pd.DataFrame()
    voting['NetMHC'] = annotation_df['NetMHC MT IC50 Score'] < 500
    voting['NetMHCpan'] = annotation_df['NetMHCpan MT IC50 Score'] < 500
    voting['MHCflurry'] = annotation_df['MHCflurry MT IC50 Score'] < 500
    voting['majority'] = voting.mode(axis=1)[0]
    annotation_df_prefiltered = annotation_df[voting['majority']]
    
    # Postfilter (some neofox feature cutoffs)

    
if __name__ == "__main__":
    main(sys.argv)