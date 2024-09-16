import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn import linear_model

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from data.get_data import get_feature_data_NEPdb

features_with_meta = get_feature_data_NEPdb()

labels = np.array(features_with_meta['response'])

feature_list = ['DAI', 'IEDB_Immunogenicity', 'MixMHCpred_score', 'PRIME_score', 'Selfsimilarity_conserved_binder', 'Tcell_predictor', 'dissimilarity_score', 'hex_alignment_score', 'recognition_potential']

features_with_meta=features_with_meta.loc[:, feature_list]
features_with_meta['Selfsimilarity_conserved_binder'] = features_with_meta['Selfsimilarity_conserved_binder'].fillna(0)

features_with_meta_array = np.array(features_with_meta)

train_features, test_features, train_labels, test_labels = train_test_split(features_with_meta_array, labels, test_size = 0.2, random_state = 42)

train_labels_int = [0 if label == 'N' else 1 for label in train_labels]
test_labels_int = [0 if label == 'N' else 1 for label in test_labels]

imputer = SimpleImputer(strategy='median')
imputer.fit(train_features)
train_features = imputer.transform(train_features)
test_features = imputer.transform(test_features)

# TODO class imbalance

rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42, min_samples_leaf=20)
rf_classifier.fit(train_features, train_labels)

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background"
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
features = features_qscore.loc[:, feature_list]
features['Selfsimilarity_conserved_binder'] = features['Selfsimilarity_conserved_binder'].fillna(0)

features = np.array(features)
qscores = np.array(qscores).flatten()

imputer = SimpleImputer(strategy='median')
imputer.fit(features)
features = imputer.transform(features)

pred_proba = rf_classifier.predict_proba(features)[:, 1]

diff = [score - pred for score, pred in zip(qscores, pred_proba)]
print("Mean", np.mean(diff))
print("Std", np.std(diff))
plt.hist(diff, bins=50)
plt.xlim([-1,1])
plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/diff_qscore.png")