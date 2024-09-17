import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

from analysis.helpers.get_data import get_relevant_features_neofox
from data.get_data import get_feature_data_NEPdb

rf = RandomForestClassifier()

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"
features_qscore = None
for cohort in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, cohort)):
        continue
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            if sample.startswith("_no"):
                continue
            tumor = sample.split("-")[0]
            normal = sample.split("-")[1]
            filename = os.path.join(directory, cohort, method, sample, "rank_sum_weighted_out_25_v2.tsv")
            if not os.path.isfile(filename):
                continue
            weighted_rank_sum_out = pd.read_csv(filename, sep="\t", header=0)
            if features_qscore is None:
                features_qscore = weighted_rank_sum_out
            else:
                features_qscore = pd.concat([features_qscore, weighted_rank_sum_out])

qscores = features_qscore.loc[:, ['qscore']]
features = features_qscore.loc[:, get_relevant_features_neofox()]
features['Selfsimilarity_conserved_binder'] = features['Selfsimilarity_conserved_binder'].fillna(0)

features = np.array(features)
qscores = np.array(qscores).flatten()

imputer = SimpleImputer(strategy='median')
imputer.fit(features)
features = imputer.transform(features)

pred_proba = rf.get_classifier().predict_proba(features)[:, 1]

diff = [score - pred for score, pred in zip(qscores, pred_proba)]
print("Mean", np.mean(diff))
print("Std", np.std(diff))
plt.hist(diff, bins=50)
plt.xlim([-1,1])
plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/diff_qscore.png")