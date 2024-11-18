import os
import pandas as pd
from scipy import stats
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot

rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 20})
rcParams['axes.unicode_minus'] = False

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/analysis/images/features/neofox/distributions"

def main(argv):
    usage = "usage: python features_neofox_plots.py -f <feature> --prefilter"
    desc = "Creates plots for neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    feature_data = get_feature_data(prefilter, cohorts)
    features = get_relevant_features()
    
    fig, axs = plt.subplots(ncols=3, nrows=4, figsize=(12,16), layout="constrained")
    
    for i, f in enumerate(features):
        if f in ["NetMHC MT IC50 Score", "NetMHCpan MT IC50 Score", "MHCflurry MT IC50 Score"]:
            continue
        data = feature_data[f].dropna()
        print(f)
        
        # descriptive statistics
        print("\t mean:", data.mean())
        print("\t median:", data.median())
        print("\t mode:", data.mode())
        print("\t std:", data.std())
        print("\t skewness:", data.skew())
        print("\t kurtosis:", data.kurtosis())
        
        pvalue_normal = stats.kstest((data - data.mean()) / data.std(), 'norm').pvalue
        
        print("\t Normal Distribution:", ("YES" if pvalue_normal >0.05 else "NO"), "with pvalue:", pvalue_normal)
        print()
        
        ax = axs[int(i/3), i%3]

        plot = qqplot(data, line ='45', ax=ax, marker='o', markerfacecolor='#94B6D2', markeredgecolor='#94B6D2')
        ax.get_lines()[1].set_color("#968C8C")
        ax.get_lines()[1].set_linewidth("2")
        ax.get_lines()[1].set_linestyle("--")
        ax.set_title(f.replace("_", " ").replace("imputed ", "").replace(" conserved binder", ""))
    fig.delaxes(axs[3][2])
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'qqplots' + ('_prefilter' if prefilter else '') + '.png'))
        
    
    
if __name__ == "__main__":
    main(sys.argv)