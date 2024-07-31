import os
import pandas as pd

# for both cohorts and each for strelka und dragen
directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background/SomaticAndTreatment/dragen"

for filename in os.listdir(directory):
    if filename == "_no_quality":
        continue
    file = os.path.join(directory, filename)
    
    # minimum 3 read count (= tumor DNA VAF * Depth)
    pvac_file = os.path.join(file, "pVACseq", "MHC_Class_I", filename.split("-")[0] + ".filtered.tsv")
    if not os.path.isfile(pvac_file):
        print("No pVACseq file")
        continue
    pvacseq = pd.read_csv(pvac_file, sep="\t", header=0)
    for _, variant in pvacseq.iterrows():
        read_count = variant['Tumor DNA VAF'] * variant['Tumor DNA Depth']
        if read_count < 3:
            print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            print(filename)
            print(variant['MT Epitope Seq'])
    