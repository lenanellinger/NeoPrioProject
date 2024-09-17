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

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features/neofox/variant_types"

def main(argv):
    usage = "usage: python features_neofox_variant_types.py --prefilter"
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features_neofox(True)
    feature_data = get_feature_data(prefilter, cohorts)
    

    for i, feature in enumerate(relevant_features):
        grouped_variants = feature_data.groupby('Variant Type')
        fig = plt.figure(figsize=(12,8))
        gs = gridspec.GridSpec(len(grouped_variants.groups.keys()), 2, width_ratios=[2, 1])
        xlim = [np.nanmin(feature_data[feature['name']]), np.nanmax(feature_data[feature['name']])]
        bin_width = (xlim[1] - xlim[0]) / 50
        
        # histogram with cohorts merged hue is variant type
        axes = []
        colors = {'missense': '#1a80bb', 'FS': '#ea801c', 'inframe_ins': '#298c8c', 'inframe_del': '#a00000'}
        for i, (name, group) in enumerate(grouped_variants):
            axes.append(plt.subplot(gs[i, 0]))
            feature_df = group.loc[:, [feature['name'], "Variant Type"]].dropna()
            feature_df_hist = feature_df
            if feature_df_hist.shape[0] == 0:
                continue
            elif np.nanmax(feature_df_hist[feature['name']]) - np.nanmin(feature_df_hist[feature['name']]) < 0.5 * bin_width:
                bin_width_small = np.nanmax(feature_df_hist[feature['name']]) - np.nanmin(feature_df_hist[feature['name']])
                if bin_width_small == 0:
                    sns.histplot(feature_df_hist, kde=True, x=feature['name'], stat='density', legend=True, color=colors[name])
                else:
                    sns.histplot(feature_df_hist, kde=True, x=feature['name'], stat='density', legend=True, color=colors[name], binwidth=bin_width_small)                    
            else:                
                sns.histplot(feature_df_hist, kde=True, x=feature['name'], stat='density', legend=True, color=colors[name], binwidth=bin_width)
            axes[i].set_title(feature['name'] + " (" + name + ")")
            axes[i].set_xlabel("mean: " + str(feature_df[feature['name']].mean()) + ", var: " + str(feature_df[feature['name']].var()))
            axes[i].set_xlim(xlim)
        
        # boxplot with cohorts merged
        feature_df = feature_data.loc[:, [feature['name'], "Variant Type"]].dropna()
        
        ax_box = plt.subplot(gs[:, 1])
        sns.boxplot(data=feature_df, x=[""] * len(feature_df), y=feature['name'], hue="Variant Type", ax=ax_box, palette=colors)

        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, feature['name'] + '_variant_types_' + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=100)
    
if __name__ == "__main__":
    main(sys.argv)