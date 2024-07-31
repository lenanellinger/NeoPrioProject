import os
import pandas as pd

directory = os.fsencode(os.path.join("_old", "strelka_hla_genotyper"))
all_alleles = []
    
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".tsv"): 
        alleles_df = pd.read_csv(os.path.join("_old", "strelka_hla_genotyper", filename), sep='\t', header=0)
        hla_alleles = [row['a1'] for _, row in alleles_df.iterrows()] + [row['a2'] for _, row in alleles_df.iterrows()]
        all_alleles.extend(hla_alleles)
dict_al = dict((x, all_alleles.count(x)) for x in sorted(set(all_alleles)))
for key in dict_al:
    print(key, ":", dict_al[key])