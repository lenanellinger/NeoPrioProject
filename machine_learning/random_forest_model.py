import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn import linear_model

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from data.get_data import get_feature_data

features_with_meta = get_feature_data()

labels = np.array(features_with_meta['response'])

features_with_meta=features_with_meta.loc[:, ['patientIdentifier', 'DAI', 'IEDB_Immunogenicity', 'MixMHCpred_score', 'PRIME_score', 'Selfsimilarity_conserved_binder', 'Tcell_predictor', 'dissimilarity_score', 'hex_alignment_score', 'recognition_potential']]
features_with_meta['Selfsimilarity_conserved_binder'] = features_with_meta['Selfsimilarity_conserved_binder'].fillna(0)

features = features_with_meta.drop(['patientIdentifier'], axis=1)
feature_list = list(features.columns)
features_with_meta_array = np.array(features_with_meta)

train_features, test_features, train_labels, test_labels = train_test_split(features_with_meta_array, labels, test_size = 0.2, random_state = 42)

train_features = train_features[:, 1:]
train_labels_int = [0 if label == 'N' else 1 for label in train_labels]
test_patients = test_features[:, 0]
test_features = test_features[:, 1:]
test_labels_int = [0 if label == 'N' else 1 for label in test_labels]

imputer = SimpleImputer(strategy='median') 
imputer.fit(train_features)
train_features = imputer.transform(train_features)
test_features = imputer.transform(test_features)

# TODO class imbalance
# TODO overfitting

rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42, min_samples_leaf=20)
rf_classifier.fit(train_features, train_labels)

clf = linear_model.Lasso(alpha=0.00001)
clf.fit(train_features, train_labels_int)

# Feature Correlation
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
corr = spearmanr(features).correlation

corr = (corr + corr.T) / 2
np.fill_diagonal(corr, 1)

# We convert the correlation matrix to a distance matrix before performing
# hierarchical clustering using Ward's linkage.
distance_matrix = 1 - np.abs(corr)
distance_matrix = np.nan_to_num(distance_matrix, nan=1)
dist_linkage = hierarchy.ward(squareform(distance_matrix))
dendro = hierarchy.dendrogram(
    dist_linkage, labels=features.columns.to_list(), ax=ax1, leaf_rotation=90
)
dendro_idx = np.arange(0, len(dendro["ivl"]))

ax2.imshow(corr[dendro["leaves"], :][:, dendro["leaves"]])
ax2.set_xticks(dendro_idx)
ax2.set_yticks(dendro_idx)
ax2.set_xticklabels(dendro["ivl"], rotation="vertical")
ax2.set_yticklabels(dendro["ivl"])
_ = fig.tight_layout()
     
plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/feature_correlation.png")  

# Feature Importances
result = permutation_importance(
    rf_classifier, test_features, test_labels, n_repeats=10, random_state=42, n_jobs=2
)

for i in result.importances_mean.argsort()[::-1]:
    if result.importances_mean[i] - 2 * result.importances_std[i] > 0:
        print(f"{feature_list[i]:<8}"
              f"{result.importances_mean[i]:.3f}"
              f" +/- {result.importances_std[i]:.3f}")
              
importances_mean = result.importances_mean
importances_std = result.importances_std

# Adjust the importance scores by dividing mean by std (to penalize high variability)
adjusted_importances = importances_mean / (importances_std + 1e-10)  # Adding a small value to avoid division by zero

# Normalize the adjusted importances to get final weights
feature_weights_adjusted = adjusted_importances / adjusted_importances.sum()

# Display the adjusted feature weights
for feature, weight in zip(features, feature_weights_adjusted):
    print(f"{feature}: Adjusted Weight = {weight:.4f}")

sorted_importances_idx = result.importances_mean.argsort()
importances = pd.DataFrame(
    result.importances[sorted_importances_idx].T,
    columns=features.columns[sorted_importances_idx],
)
ax = importances.plot.box(vert=False, whis=10)
ax.set_title("Permutation Importances (test set)")
ax.axvline(x=0, color="k", linestyle="--")
ax.set_xlabel("Decrease in accuracy score")
ax.figure.tight_layout()

plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/rf_feature_importances_test.png")

result = permutation_importance(
    rf_classifier, train_features, train_labels, n_repeats=10, random_state=42, n_jobs=2
)

sorted_importances_idx = result.importances_mean.argsort()
importances = pd.DataFrame(
    result.importances[sorted_importances_idx].T,
    columns=features.columns[sorted_importances_idx],
)
ax = importances.plot.box(vert=False, whis=10)
ax.set_title("Permutation Importances (train set)")
ax.axvline(x=0, color="k", linestyle="--")
ax.set_xlabel("Decrease in accuracy score")
ax.figure.tight_layout()

plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/rf_feature_importances_train.png")

result = permutation_importance(
    clf, test_features, test_labels_int, n_repeats=10, random_state=42, n_jobs=2
)

sorted_importances_idx = result.importances_mean.argsort()
importances = pd.DataFrame(
    result.importances[sorted_importances_idx].T,
    columns=features.columns[sorted_importances_idx],
)
ax = importances.plot.box(vert=False, whis=10)
ax.set_title("Permutation Importances (test set)")
ax.axvline(x=0, color="k", linestyle="--")
ax.set_xlabel("Decrease in accuracy score")
ax.figure.tight_layout()

plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/lasso_feature_importances_test.png")

result = permutation_importance(
    clf, train_features, train_labels_int, n_repeats=10, random_state=42, n_jobs=2
)

sorted_importances_idx = result.importances_mean.argsort()
importances = pd.DataFrame(
    result.importances[sorted_importances_idx].T,
    columns=features.columns[sorted_importances_idx],
)
ax = importances.plot.box(vert=False, whis=10)
ax.set_title("Permutation Importances (train set)")
ax.axvline(x=0, color="k", linestyle="--")
ax.set_xlabel("Decrease in accuracy score")
ax.figure.tight_layout()

plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/machine_learning/lasso_feature_importances_train.png")


# Test Set
pred_proba = rf_classifier.predict_proba(test_features)[:, 1]
pred_labels = (pred_proba > 0.5).astype(int)

precision = precision_score(test_labels_int, pred_labels)
recall = recall_score(test_labels_int, pred_labels)
f1 = f1_score(test_labels_int, pred_labels)  
auc = roc_auc_score(test_labels_int, pred_proba)

print(f"RF train accuracy: {rf_classifier.score(train_features, train_labels):.3f}")
print(f"RF test accuracy: {rf_classifier.score(test_features, test_labels):.3f}")
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc) 

pred_proba = clf.predict(test_features)
pred_labels = (pred_proba > 0.5).astype(int)

precision = precision_score(test_labels_int, pred_labels)
recall = recall_score(test_labels_int, pred_labels)
f1 = f1_score(test_labels_int, pred_labels)  
auc = roc_auc_score(test_labels_int, pred_proba)

print(f"LASSO train accuracy: {clf.score(train_features, train_labels_int):.3f}")
print(f"LASSO test accuracy: {clf.score(test_features, test_labels_int):.3f}")
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc) 



df = pd.DataFrame({'patientIdentifier': test_patients, 'pred_proba': pred_proba})
df.to_csv('/mnt/storage2/users/ahnelll1/master_thesis/training_data/NEPdb_neofox_annotations_random_forest_out.tsv', sep='\t', index=False, header=True)