from scipy.stats import ks_2samp, mannwhitneyu
import os
import pandas as pd
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/cohort_analysis/images"

def main(argv):
    usage = "usage: python cohorts_comparison_statistical_tests.py --prefilter"
    desc = "Calculaes distribution similarity tests for pVACseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens should be filtered by prefilter of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    feature_data = get_feature_data(prefilter, cohorts)
    relevant_features = get_relevant_features()
        
    f = open(os.path.join(output_dir, "cohorts_distribution_tests" + ('_prefilter' if prefilter else '') + ".txt"), "w")
    for feat in relevant_features:
        pvalue_ks = ks_2samp(feature_data[feature_data['cohort'] == 'AxelMelanomaPhD'][feat].dropna(), feature_data[feature_data['cohort'] == 'SomaticAndTreatment'][feat].dropna()).pvalue 
        pvalue_mw = mannwhitneyu(feature_data[feature_data['cohort'] == 'AxelMelanomaPhD'][feat].dropna(), feature_data[feature_data['cohort'] == 'SomaticAndTreatment'][feat].dropna()).pvalue
        f.write(feat + " KS " + str(pvalue_ks) + " " + str(pvalue_ks >= 0.05) + "\n")
        f.write(feat + " MW " + str(pvalue_mw) + " " + str(pvalue_mw >= 0.05) + "\n")
        f.write("\n")
    f.close()
     
if __name__ == "__main__":
    main(sys.argv)