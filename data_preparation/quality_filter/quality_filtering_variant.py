import os
import pandas as pd

# define directory name specific for your interest
directory = "/mnt/storage2/users/ahnelll1/master_thesis/output/SomaticAndTreatment/dragen"

"""
Checks if the variant occurs on at least 3 reads. Other variant specific quality measures are performed by pVACseq
"""
for filename in os.listdir(directory):
    if filename == "_no_quality":
        continue
    file = os.path.join(directory, filename)
    
    # minimum 3 read count (= tumor DNA VAF * Depth)
    pvac_file = os.path.join(file, "pVACseq", "MHC_Class_I", filename.split("-")[0] + ".filtered.tsv")
    if not os.path.isfile(pvac_file):
        print(f"No pVACseq file for {filename}")
        continue
    pvacseq = pd.read_csv(pvac_file, sep="\t", header=0)
    for _, variant in pvacseq.iterrows():
        read_count = variant['Tumor DNA VAF'] * variant['Tumor DNA Depth']
        if read_count < 3:
            print(f"Read Count not sufficient for sample {filename} in variant {variant['MT Epitope Seq']}")
    