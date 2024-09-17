import os
import pandas as pd

directory_ml = "/mnt/storage2/users/ahnelll1/master_thesis/output_training_data"

def get_feature_data_NEPdb():
    """
    returns a DataFrame with neofox features
    """

    features = pd.read_csv(os.path.join(directory_ml,  "NEPdb_neofox_annotations.tsv"), sep="\t", header=0)
            
    return features
