import os
import pandas as pd
from scipy import stats
from helpers.get_data import get_feature_data
import sys
from optparse import OptionParser
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features/neofox/distributions"

def main(argv):
    usage = "usage: python features_neofox_plots.py -f <feature> --prefilter"
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    parser.add_option("-f", "--feature", action="store", dest="feature", default=False, help="Feature Name")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    feature_data = get_feature_data(prefilter, cohorts)[options.feature].dropna()
   
    print((feature_data == 0).sum() / feature_data.size)
    
    # descriptive statistics
    print("mean:", feature_data.mean())
    print("median:", feature_data.median())
    print("mode:", feature_data.mode())
    print("std:", feature_data.std())
    print("skewness:", feature_data.skew())
    print("kurtosis:", feature_data.kurtosis())
    
    pvalue_normal = stats.kstest((feature_data - feature_data.mean()) / feature_data.std(), 'norm').pvalue
    
    print("Normal Distribution")
    print(("YES" if pvalue_normal >0.05 else "NO"), pvalue_normal)
    
    fig = plt.figure(figsize=(12,8))
    qqplot(feature_data, line='45')
    plt.savefig(os.path.join(output_dir, options.feature + '_qqplot' + ('_prefilter' if prefilter else '') + '.png'))
    
    
if __name__ == "__main__":
    main(sys.argv)