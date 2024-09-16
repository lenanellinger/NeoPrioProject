import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

response = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/training_data/NEPdb_neofox_annotations.tsv", sep="\t", header=0).loc[:, ["patientIdentifier", "response"]].sort_values(by=['patientIdentifier'])
rank_sum = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/training_data/NEPdb_neofox_annotations_rank_sum_weighted_out_25.tsv", sep="\t", header=0).sort_values(by=['patientIdentifier'])

# Test Set
pred_labels = [0 if qscore < 0.85 else 1 for qscore in rank_sum["qscore"]]
test_labels = [0 if label == 'N' else 1 for label in response["response"]]

accuracy = accuracy_score(test_labels, pred_labels)
precision = precision_score(test_labels, pred_labels)
recall = recall_score(test_labels, pred_labels)
f1 = f1_score(test_labels, pred_labels)
auc = roc_auc_score(test_labels, rank_sum["qscore"].to_numpy())

print("Accuracy:", accuracy)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc)