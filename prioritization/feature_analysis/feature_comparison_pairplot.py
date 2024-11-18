import os
import matplotlib.pyplot as plt
import seaborn as sns
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/feature_analysis/images"

def main(argv):
    usage = "usage: python feature_comparison_pairplot.py --prefilter"
    desc = "Creates pairplot for pVACseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens should be filtered by prefilter of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter=options.prefilter
    cohorts=['AxelMelanomaPhD', 'SomaticAndTreatment']

    features_data = get_feature_data(prefilter, cohorts)

    sns.pairplot(features_data.loc[:, get_relevant_features() + ['cohort']], hue='cohort')   
    plt.savefig(os.path.join(output_dir, 'features_pairplot' + ('_prefilter' if prefilter else '') + '.png'))
    
if __name__ == "__main__":
    main(sys.argv)