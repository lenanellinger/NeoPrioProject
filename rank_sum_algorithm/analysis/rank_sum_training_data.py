import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt

step_size = 25

response = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_neofox_annotations.tsv", sep="\t", header=0).loc[:, ["patientIdentifier", "response"]].sort_values(by=['patientIdentifier'])
rank_sum = pd.read_csv(f"/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_neofox_annotations_rank_sum_weighted_out_{step_size}.tsv", sep="\t", header=0).sort_values(by=['patientIdentifier'])

# Test Set
pred_labels = [0 if qscore < 0.5 else 1 for qscore in rank_sum["qscore"]]
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


# Calculate ROC curve
fpr, tpr, thresholds = roc_curve(test_labels, rank_sum["qscore"].to_numpy())
# Plot the ROC curve
plt.figure()
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % auc)
plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/ROC_NEPdb_{step_size}.png")