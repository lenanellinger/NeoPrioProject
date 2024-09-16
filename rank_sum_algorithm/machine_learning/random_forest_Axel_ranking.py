import pandas as pd
import numpy as np
import sys
import os

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.impute import SimpleImputer

sys.path.append(sys.path[0] + '/..')
from analysis.helpers.get_data import get_feature_data


metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
response_map = {
    'PD': 0,
    'SD': 0,
    'MR': 0,
    'CR': 1,
    'PR': 1
}

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background/AxelMelanomaPhD"
features = []
for method in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, method)):
        continue
    for sample in os.listdir(os.path.join(directory, method)):
        if sample.startswith("_no"):
            continue
        tumor = sample.split("-")[0]
        normal = sample.split("-")[1]
        filename = os.path.join(directory, method, sample, "rank_sum_weighted_out_25.tsv")
        response = metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE']
        if response.isnull().values.any() or not os.path.isfile(filename):
            continue
        weighted_rank_sum_out = pd.read_csv(filename, sep="\t", header=0)
        weighted_rank_sum_out['RESPONSE'] = response_map[response.values[0]]
        features.append(weighted_rank_sum_out[weighted_rank_sum_out['qscore'] == weighted_rank_sum_out.max()['qscore']])
features = pd.concat(features)

features_with_meta = features.loc[:,
                     ['DAI', 'IEDB_Immunogenicity', 'MixMHCpred_score', 'PRIME_score',
                      'Selfsimilarity_conserved_binder', 'Tcell_predictor', 'dissimilarity_score',
                      'hex_alignment_score', 'recognition_potential', 'imputedGeneExpression',
                      'Priority_score_imputed_fromDNA', 'Priority_score_imputed_fromRNA']]
features_with_meta['Selfsimilarity_conserved_binder'] = features_with_meta['Selfsimilarity_conserved_binder'].fillna(0)

features_with_meta_array = np.array(features_with_meta)
labels = features['RESPONSE']

train_features, test_features, train_labels, test_labels = train_test_split(features_with_meta_array, labels,
                                                                            test_size=0.2, random_state=42)

imputer = SimpleImputer(strategy='median')
imputer.fit(train_features)
train_features = imputer.transform(train_features)
test_features = imputer.transform(test_features)

# TODO class imbalance

rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42, min_samples_leaf=20)
rf_classifier.fit(train_features, train_labels)

# Test Set
pred_proba = rf_classifier.predict_proba(test_features)[:, 1]
pred_labels = (pred_proba > 0.5).astype(int)

precision = precision_score(test_labels, pred_labels)
recall = recall_score(test_labels, pred_labels)
f1 = f1_score(test_labels, pred_labels)
auc = roc_auc_score(test_labels, pred_proba)

print(f"RF train accuracy: {rf_classifier.score(train_features, train_labels):.3f}")
print(f"RF test accuracy: {rf_classifier.score(test_features, test_labels):.3f}")
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc)

"""
RF train accuracy: 0.710
RF test accuracy: 0.557
Precision: 0.375
Recall: 0.10344827586206896
F1-score: 0.16216216216216217
AUC: 0.4499579478553406
"""