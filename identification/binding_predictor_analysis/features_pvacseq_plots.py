import os
import sys
from optparse import OptionParser
from matplotlib import rc
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import seaborn as sns

from helpers.get_data import get_feature_data, get_relevant_features_pvacseq

rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 25})

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/identification/binding_predictor_analysis/images"

def main(argv):
    usage = "usage: python features_pvacseq_plots.py --prefilter"
    desc = "Creates plots for pVACseq features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens should be filtered by prefilter of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features_pvacseq()
    feature_data = get_feature_data(prefilter, cohorts)
    colors = ["#DD8047", "#A5AB81", "#94B6D2"]
    
    # IC50 histogram    
    plt.figure(figsize=(16,8))

    for feature, c in zip(relevant_features, colors): 
        feature_df = feature_data[feature].dropna()
        feature_df = feature_df[feature_df <= 1000]
        sns.kdeplot(feature_df, label=feature.replace("MT IC50 Score", ""), color=c, linewidth=4)
    plt.legend()
    plt.xlim([0, 1000])
    plt.xlabel("Mutant IC50 Score")
    plt.savefig(os.path.join(output_dir, 'binding_ic50_plots_' + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), dpi=300)
    plt.figure().clear()
    
    # venn diagram of epitopes per tool
    sets = []
    for feature in relevant_features: 
        sets.append(set(list(feature_data[feature_data[feature] < 500].index.values)))
 
    venn3(sets, tuple([f.replace("MT IC50 Score", "") for f in relevant_features]), set_colors=tuple(colors),alpha=0.7)
    plt.savefig(os.path.join(output_dir, 'binding_venn_plots_' + '_'.join(cohorts) + ('_prefilter' if prefilter else '') + '.png'), dpi=300)
    
if __name__ == "__main__":
    main(sys.argv)