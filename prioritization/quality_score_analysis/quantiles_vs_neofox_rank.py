import sys

sys.path.append(sys.path[0] + '/../..')
from helpers.get_data import get_feature_data

import json
import bisect
import numpy as np

from sklearn.metrics import mean_squared_error

from scipy.stats import pearsonr, spearmanr

step_size = 10

data = get_feature_data()
f = open(f'/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/data/quantiles_{step_size}.json')
features = json.load(f)
f.close()

for f in ['MixMHCpred', 'PRIME']:
    quants = np.array([bisect.bisect_left(features[f + "_score"]['quantiles'], num) / 100 for _, num in data[f + "_score"].items()])
    neofox_ranks = np.array([num for _, num in data[f + "_rank"].items()])
    neofox_ranks_normalized = np.array([num / max(neofox_ranks) for num in neofox_ranks])

    mask = ~np.isnan(quants) & ~np.isnan(neofox_ranks_normalized)

    # Correlation
    spearman_corr = spearmanr(quants[mask], neofox_ranks_normalized[mask])
    print("Spearman Correlation", spearman_corr)

    # Mean Square Error
    mse = mean_squared_error(quants[mask], neofox_ranks_normalized[mask])
    print("Mean Square Error", mse)