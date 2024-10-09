from scipy.stats import ks_2samp, mannwhitneyu
import os
import pandas as pd
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/analysis/images/cohorts"

def main(argv):    
    metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
        
    print(metadata["COHORT"].value_counts())
     
if __name__ == "__main__":
    main(sys.argv)