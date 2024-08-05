import os
import matplotlib.pyplot as plt
import seaborn as sns
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features"

def main(argv):
    usage = "usage: python feature_comparison_pairplot.py --majority-vote"
    desc = "Creates pairplot for pVACseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--majority-vote", action="store_true", dest="majority_vote", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    majority_vote=options.majority_vote
    cohorts=['AxelMelanomaPhD', 'SomaticAndTreatment']

    features_data = get_feature_data(majority_vote, cohorts)

    sns.pairplot(features_data.loc[:, get_relevant_features() + ['cohort']], hue='cohort')   
    plt.savefig(os.path.join(output_dir, 'features_pairplot' + ('_majority_vote' if majority_vote else '') + '.png'))
    
if __name__ == "__main__":
    main(sys.argv)