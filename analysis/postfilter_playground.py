import os
import pandas as pd
import numpy as np
from helpers.get_data import get_feature_data, get_relevant_features_neofox
import sys
from optparse import OptionParser

def main(argv):
    usage = "usage: python postfilter_playground.py --prefilter"
    desc = "Playground for postfilter cutoffs."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--prefilter", action="store_true", dest="prefilter", default=False, help="If neoantigens shouldd be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()
    
    prefilter = options.prefilter
    cohorts = ['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    relevant_features = get_relevant_features_neofox(True)
    feature_data = get_feature_data(prefilter, cohorts)
    
    print(feature_data.shape[0])
    for feature in relevant_features:
        if 'cutoff' not in feature or not feature['cutoff']['good_distribution']: # bad distribution included?
            continue
        if feature['cutoff']['use_percentage']:
            q = 0.75 # change percentage
            if feature['quantile'] == 'upper':
                q = 1-q
            cutoff = np.quantile(feature_data[feature['name']].dropna(), q)
        else:
            cutoff = feature['cutoff']['specific_value']
        if feature['quantile'] == 'lower':
            feature_data = feature_data[feature_data[feature['name']] < cutoff]
        else:
            feature_data = feature_data[feature_data[feature['name']] > cutoff]
        print(feature, cutoff)
        print(feature_data.shape[0])
        
    
    
if __name__ == "__main__":
    main(sys.argv)