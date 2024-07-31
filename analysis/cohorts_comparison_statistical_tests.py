from scipy.stats import ks_2samp, mannwhitneyu
import os
import pandas as pd
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/analysis/images/cohorts"

def main(argv):
    usage = "usage: python cohorts_comparison_statistical_tests.py --majority-vote"
    desc = "Calculaes distribution similarity tests for pVACseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--majority-vote", action="store_true", dest="majority_vote", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    majority_vote = options.majority_vote
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    feature_data = get_feature_data(majority_vote, cohorts)
    relevant_features = get_relevant_features()
        
    f = open(os.path.join(output_dir, "cohorts_distribution_tests" + ('_majority_vote' if majority_vote else '') + ".txt"), "w")
    for feat in relevant_features:
        pvalue_ks = ks_2samp(feature_data[feature_data['cohort'] == 'AxelMelanomaPhD'][feat].dropna(), feature_data[feature_data['cohort'] == 'SomaticAndTreatment'][feat].dropna()).pvalue 
        pvalue_mw = mannwhitneyu(feature_data[feature_data['cohort'] == 'AxelMelanomaPhD'][feat].dropna(), feature_data[feature_data['cohort'] == 'SomaticAndTreatment'][feat].dropna()).pvalue
        f.write(feat + " KS " + str(pvalue_ks) + " " + str(pvalue_ks >= 0.05) + "\n")
        f.write(feat + " MW " + str(pvalue_mw) + " " + str(pvalue_mw >= 0.05) + "\n")
        f.write("\n")
    f.close()
     
if __name__ == "__main__":
    main(sys.argv)