import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.impute import SimpleImputer

from data.get_data import get_feature_data

features = get_feature_data()

labels = np.array(features['response'])

features=features.loc[:, ['patientIdentifier', 'imputedGeneExpression', 'DAI', 'IEDB_Immunogenicity', 'MixMHCpred_score', 'PRIME_score', 'Selfsimilarity_conserved_binder', 'Tcell_predictor', 'dissimilarity_score', 'hex_alignment_score', 'recognition_potential']]
features['Selfsimilarity_conserved_binder'] = features['Selfsimilarity_conserved_binder'].fillna(0)

feature_list = list(features.columns)
features_array = np.array(features)

train_features, test_features, train_labels, test_labels = train_test_split(features_array, labels, test_size = 0.2, random_state = 42)

imputer = SimpleImputer(strategy='median') 
imputer.fit(train_features)
train_features = imputer.transform(train_features)
test_features = imputer.transform(test_features)

# TODO class imbalance

train_features = train_features[:, 1:]

rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
rf_classifier.fit(train_features, train_labels)

# Test Set
test_patients = test_features[:, 0]
test_features = test_features[:, 1:]
pred_proba = rf_classifier.predict_proba(test_features)[:, 1]
pred_labels = (pred_proba > 0.5).astype(int)
test_labels = [0 if label == 'N' else 1 for label in test_labels]

accuracy = accuracy_score(test_labels, pred_labels)
precision = precision_score(test_labels, pred_labels)
recall = recall_score(test_labels, pred_labels)
f1 = f1_score(test_labels, pred_labels)  
auc = roc_auc_score(test_labels, pred_proba)

print("Accuracy:", accuracy)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc) 


df = pd.DataFrame({'patientIdentifier': test_patients, 'pred_proba': pred_proba})
df.to_csv('/mnt/storage2/users/ahnelll1/master_thesis/training_data/NEPdb_neofox_annotations_random_forest_out.tsv', sep='\t', index=False, header=True)