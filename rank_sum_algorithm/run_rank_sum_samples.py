import os
from neoantigen_prioritization_rank_sum import rank_sum_qscore

"""
Runs the rank sum algorithm for all samples in the output folder
"""

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

for cohort in os.listdir(directory):
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            if sample.startswith("_no"):
                continue
            pvacseq_filename = os.path.join(directory, cohort, method, sample, "pVACseq", "MHC_Class_I", sample.split("-")[0] + ".filtered.tsv")
            neofox_filename = os.path.join(directory, cohort, method, sample, "neofox", sample.split("-")[0] + "_neofox_annotations.tsv")
            if not os.path.isfile(neofox_filename):
                continue
            
            rank_sum_qscore(pvacseq_filename, neofox_filename)
