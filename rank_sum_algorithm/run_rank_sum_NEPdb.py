import os
from neoantigen_prioritization_rank_sum import rank_sum_training_data

directory = "/mnt/storage2/users/ahnelll1/master_thesis/training_data"

neofox_filename = os.path.join(directory, "NEPdb_neofox_annotations.tsv")

rank_sum_training_data(neofox_filename)