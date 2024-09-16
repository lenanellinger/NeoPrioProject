import pandas as pd
import numpy as np
import sys

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.impute import SimpleImputer

sys.path.append(sys.path[0] + '/..')
from analysis.helpers.get_data import get_feature_data

features_with_meta_org = get_feature_data(True, ['AxelMelanomaPhD'])

metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
response_map = {
    'PD': 0,
    'SD': 0,
    'MR': 0,
    'CR': 1,
    'PR': 1
}

features_with_meta = features_with_meta_org.loc[:,
                     ['sample_id', 'DAI', 'IEDB_Immunogenicity', 'MixMHCpred_score', 'PRIME_score',
                      'Selfsimilarity_conserved_binder', 'Tcell_predictor', 'dissimilarity_score',
                      'hex_alignment_score', 'recognition_potential', 'imputedGeneExpression',
                      'Priority_score_imputed_fromDNA', 'Priority_score_imputed_fromRNA']]
features_with_meta['Selfsimilarity_conserved_binder'] = features_with_meta['Selfsimilarity_conserved_binder'].fillna(0)
features_with_meta['Priority_score_imputed_fromDNA'] = features_with_meta['Priority_score_imputed_fromDNA'].fillna(0)
features_with_meta['Priority_score_imputed_fromRNA'] = features_with_meta['Priority_score_imputed_fromRNA'].fillna(0)

aggregated_df = features_with_meta.groupby('sample_id').agg(['mean', 'median', 'std', 'max', 'min']).reset_index()
features = aggregated_df.drop(['sample_id'], axis=1)
features_with_meta_array = np.array(aggregated_df)
labels = []
for sid in aggregated_df['sample_id']:
    sample = features_with_meta_org[features_with_meta_org['sample_id'] == sid]
    if metadata[(metadata['TUMOR'] == sample['TUMOR'].values[0]) & (metadata['NORMAL'] == sample['NORMAL'].values[0])]['RESPONSE'].isnull().values.any():
        response = 0
    else:
        response = response_map[metadata[(metadata['TUMOR'] == sample['TUMOR'].values[0]) & (metadata['NORMAL'] == sample['NORMAL'].values[0])]['RESPONSE'].values[0]]
    labels.append(response)

train_features, test_features, train_labels, test_labels = train_test_split(features_with_meta_array, labels,
                                                                            test_size=0.2, random_state=42)

train_features = train_features[:, 1:]
test_patients = test_features[:, 0]
test_features = test_features[:, 1:]

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