import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats
import seaborn as sns
from helpers.get_data import get_feature_data, get_relevant_features_neofox
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images"

def main(argv):
    usage = "usage: python neoantigen_load.py --prefilter"
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens should be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features_neofox(True)
    feature_data = get_feature_data(prefilter, cohorts)
    
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,8))
    
    # after majority vote prefilter
    sizes = []
    for i in range(np.nanmax(feature_data['sample_id'])):
        if (feature_data['sample_id'] == i).sum() <= 500:
            sizes.append((feature_data['sample_id'] == i).sum())
    sizes.sort()
    sns.histplot(sizes, kde=True, bins=50, ax=axes[0])
    axes[0].set_title("after binding filter majority vote")
    axes[0].set_xlim([0, 500])
    
    # after features postfilter
    for feature in relevant_features:
        if 'cutoff' not in feature:
            continue
        if feature['quantile'] == 'lower':
            feature_data = feature_data[feature_data[feature['name']] < feature['cutoff']]
        else:
            feature_data = feature_data[feature_data[feature['name']] > feature['cutoff']]
    sizes = []
    for i in range(np.nanmax(feature_data['sample_id'])):
        if (feature_data['sample_id'] == i).sum() <= 500:
            sizes.append((feature_data['sample_id'] == i).sum())
    sizes.sort()
    sns.histplot(sizes, kde=True, bins=50, ax=axes[1])
    axes[1].set_title("after features filter")
    axes[1].set_xlim([0, 500])
    
        
    plt.savefig(os.path.join(output_dir, 'neoantigen_load_' + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=100)
    
if __name__ == "__main__":
    main(sys.argv)