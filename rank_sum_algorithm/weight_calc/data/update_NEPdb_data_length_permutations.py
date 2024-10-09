import pandas as pd
import os
from shutil import which
import matplotlib.pyplot as plt
import numpy as np

from neofox.references.references import ReferenceFolder
from neofox.model.factories import PatientFactory
from neofox.model.factories import NeoepitopeFactory
from neofox.neofox_epitope import NeoFoxEpitope
from neofox.model.conversion import ModelConverter

def get_mut_wt_peptides(pep2aln):
    # ensures that the peptides correspond even if we deal with frame shift mutations
    algn_pep_mut = pep2aln.split('|')[0].replace('\t', '').replace('"', '').replace(' ', '')
    algn_pep_wt = pep2aln.split('|')[-1].replace('\t', '').replace('"', '').replace(' ', '')
    
    result_mut = []
    result_wt = []
    for char_mut, char_wt in zip(algn_pep_mut, algn_pep_wt):
        if char_mut != '-' and char_wt != '-':  
            result_mut.append(char_mut)   
            result_wt.append(char_wt)
    return ''.join(result_mut), ''.join(result_wt)


def main():

    samples = pd.read_csv('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/data/NECID_Query.csv')
    
    print("Smaller than 8: ", (samples['antigen_len'] < 8).sum())
    print("Larger than 14: ", (samples['antigen_len'] > 14).sum())
    
    plt.hist(samples['antigen_len'], bins=40)
    plt.title("Antigen length of NEPdb dataset")
    plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/NEPdb_antigen_len.png")
    
    # filter length and MHC Class I and multiple identical samples (with different TCR information)
    samples = samples[(~pd.isna(samples['alleleA'])) & (~samples["alleleA"].str.startswith("DQA", na=True)) & (~samples["alleleA"].str.startswith("DPA", na=True))]
    samples = samples.drop_duplicates(subset=['mut_peptide', 'wt_peptide', 'genesymbol', 'alleleA', 'Tumor Type'], keep='first')
    
    samples['custom_id'] = samples.groupby(['mut_peptide', 'wt_peptide']).ngroup() + 1
    
    # if length is greater than 14, we extract all 9mers around the variation and calculate features of each 9mer. For training only the best 9mer is kept
    samples_updated = []
    for ind, row in samples.iterrows():
        print(ind)
        if row['antigen_len'] > 14:
            if pd.isna(row['pep2aln']):
                continue
                
            pep_mut, pep_wt = get_mut_wt_peptides(row['pep2aln'])
            if len(pep_mut) != len(pep_wt):
                raise ValueError()
            alignment_str = row['pep2aln'][row['pep2aln'].find('|'):row['pep2aln'].rfind('|')+1]

            pos_first_mut = alignment_str.find('.')
            pos_last_mut = alignment_str.rfind('.')

            if pos_last_mut-pos_first_mut >= 9:
                row['mut_peptide'] = pep_mut[pos_first_mut:pos_last_mut+1]
                row['wt_peptide'] = pep_wt[pos_first_mut:pos_last_mut+1]
                row['antigen_len'] = pos_last_mut - pos_first_mut
                samples_updated.append(row)
            else:
                substrings_mut = [pep_mut[i:i+9] for i in range(max(0, pos_last_mut - 8), min(pos_first_mut + 1, len(pep_mut) - 9 + 1))]
                substrings_wt = [pep_wt[i:i+9] for i in range(max(0, pos_last_mut - 8), min(pos_first_mut + 1, len(pep_wt) - 9 + 1))]
                for s_mut, s_wt in zip(substrings_mut, substrings_wt):
                    row_copy = row.copy() 
                    row_copy['mut_peptide'] = s_mut
                    row_copy['wt_peptide'] = s_wt
                    row_copy['antigen_len'] = 9
                    samples_updated.append(row_copy)
        else:
            samples_updated.append(row)
    
    
    samples_updated_df = pd.DataFrame(samples_updated)
    samples_updated_df.to_csv('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/data/NECID_Query_updated.csv', index=False)
    

if __name__ == "__main__":
    main()