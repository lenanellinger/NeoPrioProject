import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, confusion_matrix
from sklearn.calibration import calibration_curve
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np
from scipy.stats import pearsonr, spearmanr
import seaborn as sns

rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 20})
rcParams['axes.unicode_minus'] = False

step_size = 10

response = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox.tsv", sep="\t", header=0).loc[:, ["patientIdentifier", "response"]]
prediction = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox_xgb_out.tsv", sep="\t", header=0)
rank_sum = pd.read_csv(f"/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox_rank_sum_weighted_out_{step_size}.tsv", sep="\t", header=0).loc[:, ["patientIdentifier", "qscore"]]

comparison_df = prediction.set_index('patientIdentifier').join(rank_sum.set_index('patientIdentifier'), how='inner').join(response.set_index('patientIdentifier'), how='inner')
comparison_df["response"] = [0 if r == 'N' else 1 for r in comparison_df["response"]]

# rank sum qscore vs true label

# find threshold on train set
tmp_df = rank_sum.set_index('patientIdentifier').join(response.set_index('patientIdentifier'), how='inner')
train_set = tmp_df[~tmp_df.index.isin(comparison_df.index)]
train_set["response"] = [0 if r == 'N' else 1 for r in train_set["response"]]

# find threshold of qscore by maximizing f1
best_f1 = 0
best_threshold = 0
f1_scores = []

for threshold in np.arange(0, 1, 0.01).tolist():
    y_pred = (train_set["qscore"].to_numpy() >= threshold).astype(int)
    f1 = f1_score(train_set["response"].to_numpy(), y_pred)
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
plt.close()

print("Best Threshold:", best_threshold)
print("Best F1 Score:", best_f1)

# Youden’s J Statistic
fpr, tpr, thresholds = roc_curve(train_set["response"], train_set["qscore"])

youdens_j = tpr - fpr
best_threshold_y = thresholds[youdens_j.argmax()]

print("Best Threshold based on Youden’s J:", best_threshold_y)

# statistics on test set
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
plt.close()

# rank sum qscore vs prediction 
# ROC curve
fpr1, tpr1, _ = roc_curve(comparison_df["response"], comparison_df["qscore"])
auc1 = roc_auc_score(comparison_df["response"], comparison_df["qscore"])

fpr2, tpr2, _ = roc_curve(comparison_df["response"], comparison_df["pred_proba"])
auc2 = roc_auc_score(comparison_df["response"], comparison_df["pred_proba"])

plt.figure(figsize=(8, 6))
plt.plot(fpr1, tpr1, label=f'Quality Score (area = {auc1:.2f})', color="#A5AB81", linewidth=5)
plt.plot(fpr2, tpr2, label=f'XGBoost Classifier (area = {auc2:.2f})', color="#94B6D2", linewidth=5)
plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/ROC_NEPdb_RS_{step_size}.png", dpi=300)
plt.close()

# calibration curve
true_pos1, pred_pos1 = calibration_curve(comparison_df["response"], comparison_df["qscore"], n_bins=10)
true_pos2, pred_pos2 = calibration_curve(comparison_df["response"], comparison_df["pred_proba"], n_bins=10)
 
plt.plot(pred_pos1, true_pos1, marker='o', linewidth=1, label='Quality Score')
plt.plot(pred_pos2, true_pos2, marker='o', linewidth=1, label='XGBoost Classifier')
plt.plot([0, 1], [0, 1], linestyle='--', label='Perfectly Calibrated')
plt.xlabel('Predicted Probability')
plt.ylabel('True Probability')
plt.legend(loc='best')
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/CALIBR_NEPdb_RS_{step_size}.png", dpi=300)
plt.close()

# Precision Recall Curve
precision1, recall1, _ = precision_recall_curve(comparison_df["response"], comparison_df["qscore"])
ap1 = average_precision_score(comparison_df["response"], comparison_df["qscore"])

precision2, recall2, _ = precision_recall_curve(comparison_df["response"], comparison_df["pred_proba"])
ap2 = average_precision_score(comparison_df["response"], comparison_df["pred_proba"])

plt.figure(figsize=(8, 6))
plt.plot(recall1, precision1, label=f'Rank Sum QScore (AP = {ap1:.2f})')
plt.plot(recall2, precision2, label=f'XGBoost Prediction Probability (AP = {ap2:.2f})')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/PR_NEPdb_RS_{step_size}.png", dpi=300)
plt.close()

# correlation
correlation = spearmanr(comparison_df['qscore'], comparison_df['pred_proba'])

print(f"Spearman Correlation: {correlation}")

# Scatter plot
plt.figure(figsize=(6, 6))
plt.scatter(comparison_df["qscore"], comparison_df["pred_proba"], alpha=0.5)
plt.plot([min(comparison_df["qscore"]), max(comparison_df["qscore"])], [min(comparison_df["qscore"]), max(comparison_df["qscore"])], color='red', linestyle='--')  # y=x line
plt.xlabel("Rank Sum Quality Score")
plt.ylabel("XGBoost Prediction Probability")
plt.title(f"Scatter Plot (Spearman r = {correlation.statistic:.2f})")
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/qscore_xgb_NEPdb_scatter_{step_size}.png")

# Bland-Altman Plot
mean_pred = (comparison_df["qscore"] + comparison_df["pred_proba"]) / 2
diff_pred = comparison_df["qscore"] - comparison_df["pred_proba"]
mean_diff = np.mean(diff_pred)
std_diff = np.std(diff_pred)

plt.figure(figsize=(6, 6))
plt.scatter(mean_pred, diff_pred, alpha=0.5)
plt.axhline(mean_diff, color='red', linestyle='--')
plt.axhline(mean_diff + 1.96 * std_diff, color='blue', linestyle='--')
plt.axhline(mean_diff - 1.96 * std_diff, color='blue', linestyle='--')
plt.xlabel("Mean of Predictions")
plt.ylabel("Difference between Predictions")
plt.title("Bland-Altman Plot")
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/qscore_xgb_NEPdb_bland_altman_{step_size}.png")


