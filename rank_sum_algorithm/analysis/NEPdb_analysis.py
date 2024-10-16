import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, confusion_matrix
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns

step_size = 25

response = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox.tsv", sep="\t", header=0).loc[:, ["patientIdentifier", "response"]]
prediction = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox_xgb_out.tsv", sep="\t", header=0)
rank_sum = pd.read_csv(f"/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox_rank_sum_weighted_out_{step_size}.tsv", sep="\t", header=0).loc[:, ["patientIdentifier", "qscore"]]

# rank sum qscore vs true label
comparison_df = response.set_index('patientIdentifier').join(rank_sum.set_index('patientIdentifier'), how='inner')
comparison_df["response"] = [0 if r == 'N' else 1 for r in comparison_df["response"]]

# find threshold of qscore by maximizing f1
best_f1 = 0
best_threshold = 0
f1_scores = []

for threshold in np.arange(0, 1, 0.01).tolist():
    y_pred = (comparison_df["qscore"].to_numpy() >= threshold).astype(int)
    f1 = f1_score(comparison_df["response"].to_numpy(), y_pred)
    f1_scores.append(f1)
    if f1 > best_f1:
        best_f1 = f1
        best_threshold = threshold

plt.figure(figsize=(10, 7))
plt.plot(np.arange(0, 1, 0.01).tolist(), f1_scores)
plt.xlabel('Threshold')
plt.ylabel('F1 Score')
plt.title('Rank Sum QScore on NEPdb dataset Classification Threshold')
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/F1_RS_{step_size}.png")

print("Best Threshold:", best_threshold)
print("Best F1 Score:", best_f1)

pred_rs_labels = (comparison_df["qscore"].to_numpy() >= best_threshold).astype(int)
true_labels = comparison_df["response"].to_numpy()

accuracy = accuracy_score(true_labels, pred_rs_labels)
precision = precision_score(true_labels, pred_rs_labels)
recall = recall_score(true_labels, pred_rs_labels)
f1 = f1_score(true_labels, pred_rs_labels)
auc = roc_auc_score(true_labels, comparison_df["qscore"].to_numpy())

print("Accuracy RS QScore:", accuracy)
print("Precision RS QScore:", precision)
print("Recall RS QScore:", recall)
print("F1-score RS QScore:", f1)
print("AUC RS QScore:", auc)

cm = confusion_matrix(true_labels, pred_rs_labels, labels=[0, 1])
plt.figure(figsize=(10, 7))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=[0, 1], yticklabels=[0, 1])
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title('Confusion Matrix for Rank Sum Quality Score Classifier')
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/CM_RS_{step_size}.png")

# rank sum qscore vs prediction 
# ROC curve
comparison_df = prediction.set_index('patientIdentifier').join(rank_sum.set_index('patientIdentifier'), how='inner').join(response.set_index('patientIdentifier'), how='inner')
comparison_df["response"] = [0 if r == 'N' else 1 for r in comparison_df["response"]]

fpr1, tpr1, _ = roc_curve(comparison_df["response"], comparison_df["qscore"])
auc1 = roc_auc_score(comparison_df["response"], comparison_df["qscore"])

fpr2, tpr2, _ = roc_curve(comparison_df["response"], comparison_df["pred_proba"])
auc2 = roc_auc_score(comparison_df["response"], comparison_df["pred_proba"])

plt.figure(figsize=(8, 6))
plt.plot(fpr1, tpr1, label=f'Rank Sum QScore (AUC = {auc1:.2f})')
plt.plot(fpr2, tpr2, label=f'XGBoost Prediction Probability (AUC = {auc2:.2f})')
plt.plot([0, 1], [0, 1], 'k--', label='Random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/ROC_NEPdb_RS_{step_size}.png", dpi=300)

# Precision Recall Curve
precision1, recall1, _ = precision_recall_curve(comparison_df["response"], comparison_df["qscore"])
ap1 = average_precision_score(comparison_df["response"], comparison_df["qscore"])

# Compute Precision-Recall curve and Average Precision for Classifier 2
precision2, recall2, _ = precision_recall_curve(comparison_df["response"], comparison_df["pred_proba"])
ap2 = average_precision_score(comparison_df["response"], comparison_df["pred_proba"])

# Plot Precision-Recall curves
plt.figure(figsize=(8, 6))
plt.plot(recall1, precision1, label=f'Rank Sum QScore (AP = {ap1:.2f})')
plt.plot(recall2, precision2, label=f'XGBoost Prediction Probability (AP = {ap2:.2f})')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/PR_NEPdb_RS_{step_size}.png", dpi=300)

# correlation
correlation = pearsonr(comparison_df['qscore'], comparison_df['pred_proba'])

print(f"Pearson Correlation: {correlation}")


