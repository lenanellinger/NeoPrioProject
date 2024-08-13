import os
import numpy as np
import json
import sys
from optparse import OptionParser

sys.path.append(sys.path[0] + '/..')
from analysis.helpers.get_data import get_relevant_features, get_feature_data


def main(argv):  
    feature_data = get_feature_data()
    
    print(feature_data.isna().sum()[feature_data.isna().sum() != 0] / feature_data.shape[0])
        

    
if __name__ == "__main__":
    main(sys.argv)