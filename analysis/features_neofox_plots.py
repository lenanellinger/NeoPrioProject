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
from matplotlib import rc, rcParams

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features/neofox"

def main(argv):
    usage = "usage: python features_neofox_plots.py --prefilter --use-interval"
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by binding prefilter")
    parser.add_option("--use-interval", action="store_true", dest="use_interval", default=False, help="If restricted intervals should be used to plot")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features_neofox(True)
    feature_data = get_feature_data(prefilter, cohorts)

    for i, feature in enumerate(relevant_features):
        fig = plt.figure(figsize=(12,8))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 2, 1])
            
        # violin plot with both cohorts
        ax_violin = plt.subplot(gs[0])
        features_df_violin = feature_data.loc[:, [feature['name'], 'cohort']]
        features_df_violin = features_df_violin.rename(columns={feature['name']: 'value'})
        features_df_violin['name'] = feature['name']
        sns.violinplot(ax=ax_violin, data=features_df_violin, x="name", y='value', hue="cohort", split=True, inner="quart")
        
        # histogram with cohorts merged
        ax_hist = plt.subplot(gs[1])
        feature_df = feature_data[feature['name']].dropna()
        if options.use_interval and 'interval' in feature:
            feature_df_hist = feature_df[(feature['interval'][0] <= feature_df) & (feature_df <= feature['interval'][1])]
        else:
            feature_df_hist = feature_df
        sns.histplot(feature_df_hist, kde=True, bins=50)
        for q, c in zip([0.01, 0.05, 0.1, 0.25, 0.5, 0.75], ['#a70000', '#ff0000', '#ff5252', '#ff7b7b', '#ffbaba', '#ffd5d5']):
            if feature['quantile'] == 'upper':
                q = 1-q
            quant = np.quantile(feature_df, q)
            plt.axvline(x = quant, color = c, label = str(q*100) + '% = ' + str(quant)) 
        if options.use_interval and 'interval' in feature:
            ax_hist.set_xlim(feature['interval'])
        ax_hist.set_title(feature['name'])
        ax_hist.set_xlabel("mean: " + str(feature_df.mean()) + ", var: " + str(feature_df.var()))
        ax_hist.legend()
        
        # boxplot with cohorts merged
        ax_box = plt.subplot(gs[2])
        ax_box.boxplot(feature_df)

        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, feature['name'] + '_plots_' + ('interval_' if options.use_interval and 'interval' in feature else '') + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=100)
        
    rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 30})
    rcParams['axes.unicode_minus'] = False
    for i, feature in enumerate(relevant_features):
        plt.figure(figsize=(12,8))
        
        # histogram 
        feature_df = feature_data[feature['name']].dropna()
        print(feature['name'])
        print(min(feature_df))
        print(max(feature_df))
        if options.use_interval and 'interval' in feature:
            feature_df_hist = feature_df[(feature['interval'][0] <= feature_df) & (feature_df <= feature['interval'][1])]
        else:
            feature_df_hist = feature_df
        sns.histplot(feature_df_hist, kde=True, bins=50, color="#94B6D2")
        for q, c in zip([0.25, 0.5, 0.75], ['#ff0000', '#ff7b7b', '#ffd5d5']):
            if feature['quantile'] == 'upper':
                q = 1-q
            quant = np.quantile(feature_df, q)
            plt.axvline(x = quant, color = c, label = str(q*100) + '% = ' + str(round(quant*100)/100), linewidth=4) 
        if options.use_interval and 'interval' in feature:
            plt.xlim(feature['interval'])
        plt.ylabel("#neoepitopes")
        plt.xlabel(feature['name'].replace("_", " "))
        plt.legend()

        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, "histograms", feature['name'] + '_histogram_' + ('interval_' if options.use_interval and 'interval' in feature else '') + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=100)
    
if __name__ == "__main__":
    main(sys.argv)