import os
import pandas as pd
import itertools
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser
import seaborn as sns
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 15})
rcParams['axes.unicode_minus'] = False
  

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features"

def main(argv):
    usage = "usage: python feature_comparison_correlations.py --prefilter"
    desc = "Calculates correlation coefficients for pVACseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features()
    feature_data = get_feature_data(prefilter, cohorts)
    rel_features_df = feature_data.loc[:, relevant_features]
    rel_features_df.columns = rel_features_df.columns.str.replace('_', ' ')
    corr = spearmanr(rel_features_df, nan_policy='omit').statistic
    names = [f.replace('_', ' ') for f in relevant_features]
    corr_df = pd.DataFrame(corr, columns=names, index=names)
    
    heatmap = sns.heatmap(corr_df, cmap=sns.diverging_palette(42, 207, 60, 60, as_cmap=True), vmin=-1, vmax=1)
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')
    
    plt.savefig(os.path.join(output_dir, 'feature_correlation' + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=300)
    plt.close()
    
    f = open(os.path.join(output_dir, "features_correlation_coefficients" + ('_prefilter' if prefilter else '') + ".txt"), "w")
    for cohort in cohorts:
        f.write(cohort + "\n")
        features = feature_data[feature_data['cohort'] == cohort]
        
        for f1, f2 in itertools.combinations(relevant_features, 2):
            corr_p = features[f1].corr(features[f2], method='pearson') 
            corr_k = features[f1].corr(features[f2], method='kendall')
            corr_s = features[f1].corr(features[f2], method='spearman') 
            if corr_p >= 0.4:
                f.write("\t" + f1 + " " + f2 + " pearson correlation: " + str(corr_p) + "\n")
            if corr_k >= 0.4:
                f.write("\t" + f1 + " " + f2 + " kendall correlation: " + str(corr_k) + "\n")
            if corr_s >= 0.4:
                f.write("\t" + f1 + " " + f2 + " spearman correlation: " + str(corr_s) + "\n")

    f.write("merged\n")        
    for f1, f2 in itertools.combinations(relevant_features, 2):
        corr_s = spearmanr(feature_data[f1], feature_data[f2], nan_policy='omit')
        if corr_s.pvalue < 0.05 and abs(corr_s.statistic) >= 0.4:
            f.write(f"\t{f1} vs {f2} Spearman Correlation (at least moderate): {corr_s.statistic}, p-value {corr_s.pvalue}\n")
    f.close()
    
if __name__ == "__main__":
    main(sys.argv)