import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot_2samples
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser
  
output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/cohort_analysis/images"

def main(argv):
    usage = "usage: python cohorts_comparison_qq_plot.py --prefilter"
    desc = "Creates pairplot for pVACseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens should be filtered by prefilter of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features()
    feature_data = get_feature_data(prefilter, cohorts)

    fig = plt.figure(figsize=(30,8))
    gs = gridspec.GridSpec(2, int(len(relevant_features) / 2) )

    for i, feature in enumerate(relevant_features):
        ax = plt.subplot(gs[i])
        
        coh_0 = feature_data[feature_data['cohort'] == cohorts[0]][feature]
        coh_1 = feature_data[feature_data['cohort'] == cohorts[1]][feature]
        
        pp_x = sm.ProbPlot(coh_0)
        pp_y = sm.ProbPlot(coh_1)
        qqplot_2samples(pp_x, pp_y, line='45', ax=ax, xlabel=(cohorts[0] if len(coh_0) > len(coh_1) else cohorts[1]), ylabel=(cohorts[0] if len(coh_0) < len(coh_1) else cohorts[1]))
        ax.set_title(feature)
        
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cohorts_qq' + ('_prefilter' if prefilter else '') + '.png'))
    plt.figure().clear()
    
if __name__ == "__main__":
    main(sys.argv)