import os
import pandas as pd

"""
Combines all pVACseq and neofox annotations per cohort into one file
"""

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

sample_id = 0
for cohort in os.listdir(directory):
    df = None
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            if sample.startswith("_no"):
                continue
            neofox_filename = os.path.join(directory, cohort, method, sample, "neofox", sample.split("-")[0] + "_neofox_annotations.tsv")
            if not os.path.isfile(neofox_filename):
                sample_id += 1
                continue
            neofox_df = pd.read_csv(neofox_filename, sep="\t", header=0)
            pvacseq_filename = os.path.join(directory, cohort, method, sample, "pVACseq", "MHC_Class_I", sample.split("-")[0] + ".filtered.tsv")
            pvacseq_df = pd.read_csv(pvacseq_filename, sep="\t", header=0)
            if neofox_df.shape[0] != pvacseq_df.shape[0]:
                print("AAAAA", cohort, method, sample)
            neofox_df_merge = neofox_df.merge(right=pvacseq_df.loc[:, ['Chromosome', 'Start', 'Stop', 'Transcript', 'Variant Type', 'MHCflurry WT IC50 Score', 'MHCflurry MT IC50 Score', 'MHCflurry WT Percentile', 'MHCflurry MT Percentile', 
            'NetMHC WT IC50 Score', 'NetMHC MT IC50 Score', 'NetMHC WT Percentile', 'NetMHC MT Percentile', 
            'NetMHCpan WT IC50 Score', 'NetMHCpan MT IC50 Score', 'NetMHCpan WT Percentile',  'NetMHCpan MT Percentile', 
            'Predicted Stability', 'Half Life', 'Stability Rank', 'NetMHCstab allele']], on=['Chromosome', 'Start', 'Stop', 'Transcript'])
            if neofox_df.shape[0] != neofox_df_merge.shape[0]:
                print(cohort, method, sample)
                print(neofox_df.loc[:, ['Chromosome', 'Start', 'Stop', 'Transcript']])
            neofox_df_merge['sample_id'] = sample_id
            neofox_df_merge['TUMOR'] = sample.split("-")[0]
            neofox_df_merge['NORMAL'] = sample.split("-")[1]
            sample_id += 1
            if df is None:
                df = neofox_df_merge
            else:
                df = pd.concat([df, neofox_df_merge])
    print(df)
    df.to_csv(os.path.join(directory, cohort, cohort + "_all.tsv"), sep="\t", index=False)
           