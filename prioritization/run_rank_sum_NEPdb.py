import os
from neoantigen_prioritization_rank_sum import rank_sum_qscore_training_data

"""
Runs the rank sum algorithm for the NEPdb data
"""

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_training_data"

neofox_filename = os.path.join(directory, "NEPdb_filtered_neofox.tsv")

rank_sum_qscore_training_data(neofox_filename, 25)