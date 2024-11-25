import pandas as pd
import os

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/all"

best_peptides = []

for filename in os.listdir(directory):
    neofox_annotation = pd.read_csv(os.path.join(directory, filename), sep="\t", header=0)
    neofox_annotation = neofox_annotation[neofox_annotation['affinityMutated'] != 0.0]
    best_peptides.append(neofox_annotation[neofox_annotation['affinityMutated']==neofox_annotation['affinityMutated'].min()].head(1))

samples = pd.read_csv('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/data/NECID_Query.csv', header=0)
print(samples.shape[0])
samples = samples[(~pd.isna(samples['alleleA'])) & (~samples["alleleA"].str.startswith("DQA", na=True)) & (~samples["alleleA"].str.startswith("DPA", na=True))]
print(samples.shape[0])
samples = samples.drop_duplicates(subset=['mut_peptide', 'wt_peptide', 'genesymbol', 'alleleA', 'Tumor Type'], keep='first')
print(samples.shape[0])
samples = samples.drop_duplicates(subset=['mut_peptide', 'wt_peptide', 'genesymbol', 'Tumor Type'], keep='first')
print(samples.shape[0])
samples = samples.drop_duplicates(subset=['mut_peptide', 'wt_peptide', 'genesymbol'], keep='first')
print(samples.shape[0])
samples = samples.drop_duplicates(subset=['mut_peptide', 'wt_peptide'], keep='first')
print(samples.shape[0])
best_peptides_df = pd.concat(best_peptides)
best_peptides_df.to_csv("/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox.tsv", sep="\t", index=False)