import os
import shutil

"""
Deletes all neofox output folders in output folder structure
"""

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

for cohort in os.listdir(directory):
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            if sample.startswith("_no") or not os.path.isdir(os.path.join(directory, cohort, method, sample, 'neofox')):
                continue
            shutil.rmtree(os.path.join(directory, cohort, method, sample, 'neofox'))