import pandas as pd
import sys

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/cohort_analysis/images"

def main(argv):    
    metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
        
    print(metadata["COHORT"].value_counts())
     
if __name__ == "__main__":
    main(sys.argv)