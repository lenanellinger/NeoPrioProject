import os
import itertools
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser
import seaborn as sns
import matplotlib.pyplot as plt
  

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
    
    sns.heatmap(feature_data.loc[:, relevant_features].corr(numeric_only=True), cmap="YlGnBu")

    plt.savefig(os.path.join(output_dir, 'feature_correlation' + ('_prefilter' if prefilter else '') + '.png'), bbox_inches='tight', dpi=100)
    
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
        corr_p = feature_data[f1].corr(feature_data[f2], method='pearson') 
        corr_k = feature_data[f1].corr(feature_data[f2], method='kendall')
        corr_s = feature_data[f1].corr(feature_data[f2], method='spearman') 
        if corr_p >= 0.4:
            f.write("\t" + f1 + " " + f2 + " pearson correlation: " + str(corr_p) + "\n")
        if corr_k >= 0.4:
            f.write("\t" + f1 + " " + f2 + " kendall correlation: " + str(corr_k) + "\n")
        if corr_s >= 0.4:
            f.write("\t" + f1 + " " + f2 + " spearman correlation: " + str(corr_s) + "\n")
    f.close()
    
if __name__ == "__main__":
    main(sys.argv)