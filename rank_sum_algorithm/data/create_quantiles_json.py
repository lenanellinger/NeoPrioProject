import numpy as np
import json
import sys
from optparse import OptionParser

sys.path.append(sys.path[0] + '/../..')
from analysis.helpers.get_data import get_relevant_features_neofox, get_feature_data


def main():
    usage = "usage: python create_quantiles_json.py --step 10"
    desc = "Creates quantile thresholds for each feature in the desired step width."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--step", action="store", dest="step", help="Step size of quantiles as percentage value")
    (options, args) = parser.parse_args()
    
    step = int(options.step)
    
    relevant_features = get_relevant_features_neofox(True)
    feature_data = get_feature_data()   
    
    result = {}
    
    for feature in relevant_features:
        result[feature['name']] = {
                'direction': feature['quantile'],
                'quantiles': []
            }
        for q in range(step, 100, step):
            quant = np.nanquantile(feature_data[feature['name']], q/100)
            result[feature['name']]['quantiles'].append(quant)
    with open('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/data/quantiles_' + options.step + '.json', 'w') as f:
        json.dump(result, f, ensure_ascii=False, indent=4)

    
if __name__ == "__main__":
    main()