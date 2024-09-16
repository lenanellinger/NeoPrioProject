import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import seaborn as sns
import sys
from optparse import OptionParser

sys.path.append(sys.path[0] + '/../..')
from analysis.helpers.get_data import get_relevant_features_neofox
from rank_sum_algorithm.machine_learning.data.get_data import get_feature_data_NEPdb

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/analysis/images/features/neofox"

def main(argv):
    usage = "usage: python features_neofox_plots.py --use-interval"
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--use-interval", action="store_true", dest="use_interval", default=False, help="If restricted intervals should be used to plot")
    (options, args) = parser.parse_args()
    
    relevant_features = get_relevant_features_neofox(True)
    feature_data = get_feature_data_NEPdb()

    for i, feature in enumerate(relevant_features):
        fig = plt.figure(figsize=(12,8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
        
        # histogram 
        ax_hist = plt.subplot(gs[0])
        if not feature['name'] in feature_data:
            continue
        feature_df = feature_data[feature['name']].dropna()
        if options.use_interval and 'interval' in feature:
            feature_df_hist = feature_df[(feature['interval'][0] <= feature_df) & (feature_df <= feature['interval'][1])]
        else:
            feature_df_hist = feature_df
        sns.histplot(feature_df_hist, kde=True, bins=50)
        for q, c in zip([0.01, 0.05, 0.1, 0.25, 0.5, 0.75], ['#a70000', '#ff0000', '#ff5252', '#ff7b7b', '#ffbaba', '#ffd5d5']):
            if feature['quantile'] == 'upper':
                q = 1-q
            if feature_df.shape[0] != 0:
                quant = np.quantile(feature_df, q)
                plt.axvline(x = quant, color = c, label = str(q*100) + '% = ' + str(quant)) 
        if options.use_interval and 'interval' in feature:
            ax_hist.set_xlim(feature['interval'])
        ax_hist.set_title(feature['name'])
        ax_hist.set_xlabel("mean: " + str(feature_df.mean()) + ", var: " + str(feature_df.var()))
        ax_hist.legend()
        
        # boxplot with cohorts merged
        ax_box = plt.subplot(gs[1])
        ax_box.boxplot(feature_df)

        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, feature['name'] + '_plots_' + ('interval_' if options.use_interval and 'interval' in feature else '') + '.png'), bbox_inches='tight', dpi=100)
    
if __name__ == "__main__":
    main(sys.argv)