import pandas as pd
import shutil
import os

metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
duplicates = metadata[(metadata["COHORT"] + metadata["PATIENT"]).duplicated('last')]

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background/AxelMelanomaPhD/strelka"
print(duplicates)
for _, d in duplicates.iterrows():
    file = os.path.join(directory, d['TUMOR'] + "-" + d['NORMAL'])
    if os.path.isdir(file):
        print("duplicated: move", file)
        shutil.move(file, os.path.join(directory, "_no_quality", d['TUMOR'] + "-" + d['NORMAL']))
    else: 
        print("duplicated but no result dir:", file)