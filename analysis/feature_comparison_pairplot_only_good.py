from matplotlib import gridspec
import numpy as np
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser
from itertools import combinations
import matplotlib.pyplot as plt
import os

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features"

def main(argv):
    usage = "usage: python feature_comparison_pairplot_only_good.py --majority-vote"
    desc = "Creates pairplot for pvacseq and neofox features which only plots sample which have a good value in both features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-p", "--percentage", action="store", dest="percentage", default=10, type=int, help="Qunaile percentage used to filter")
    parser.add_option("--majority-vote", action="store_true", dest="majority_vote", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    quantile = options.percentage / 100
    majority_vote = options.majority_vote
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features(True)
    feature_data = get_feature_data(majority_vote, cohorts)
    fig = plt.figure(figsize=(20,20))
    gs = gridspec.GridSpec(len(relevant_features)-1, len(relevant_features)-1)
    row = 0
    column = 0
    for f1, f2 in combinations(relevant_features, 2):
        ax = plt.subplot(gs[row*(len(relevant_features)-1)+column])
        q1 = np.quantile(feature_data[f1['name']], quantile if f1['quantile'] == 'lower' else 1-quantile)
        q2 = np.quantile(feature_data[f2['name']], quantile if f2['quantile'] == 'lower' else 1-quantile)
        
        good_1 = feature_data[f1['name']] < q1 if f1['quantile'] == 'lower' else feature_data[f1['name']] > q1
        good_2 = feature_data[f2['name']] < q2 if f2['quantile'] == 'lower' else feature_data[f2['name']] > q2
        good_samples = feature_data[good_1 & good_2]
    
        ax.scatter(good_samples[f2['name']], good_samples[f1['name']], s=0.5)
        if column == row:
            ax.set_xlabel(f2['name'])
            ax.set_ylabel(f1['name'])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
            
        if column == len(relevant_features)-2:
            column = row + 1
            row += 1
        else:
            column += 1
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'features_pairplot_good_' + str(options.percentage) + ('_majority_vote' if majority_vote else '') + '.png'))
    
if __name__ == "__main__":
    main(sys.argv)