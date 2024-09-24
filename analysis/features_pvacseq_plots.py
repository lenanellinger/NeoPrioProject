import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from helpers.get_data import get_feature_data, get_relevant_features_pvacseq
import sys
from optparse import OptionParser
import seaborn as sns

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/analysis/images/features/pvacseq"

def main(argv):
    usage = "usage: python features_pvacseq_plots.py --prefilter"
    desc = "Creates plots for pVACseq features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features_pvacseq()
    feature_data = get_feature_data(prefilter, cohorts)

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(12,8))
    
    for i, feature in enumerate(relevant_features): 
        feature_df = feature_data[feature].dropna()
        feature_df = feature_df[feature_df <= 1000]
        sns.histplot(feature_df, kde=True, bins=50, ax=axes[i])
        axes[i].set_title(feature)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'binding_plots_' + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=100)
    plt.figure().clear()
    
if __name__ == "__main__":
    main(sys.argv)