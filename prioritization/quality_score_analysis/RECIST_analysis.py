import os
import sys 
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from scipy.stats import wilcoxon

sys.path.append(sys.path[0] + '/../..')
from helpers.get_data import get_relevant_features_neofox
from prioritization.weight_calculation.xgboost_classifier_model import XGBoostClassifierModel

rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 20})
rcParams['axes.unicode_minus'] = False

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output/AxelMelanomaPhD"

metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
response_map = {
    'PD': 0,
    'SD': 0,
    'MR': 0,
    'CR': 1,
    'PR': 1
}

output_df = []

topx = 10
step_size = 10

rf = XGBoostClassifierModel()

def sum_top_x(df, column, number):
    sorted_scores = df[column].sort_values(ascending=False)
    top_10_scores = sorted_scores.head(number)
    
    return top_10_scores.median()

for method in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, method)):
        continue
    for sample in os.listdir(os.path.join(directory, method)):
        if sample.startswith("_no"):
            continue
        tumor = sample.split("-")[0]
        normal = sample.split("-")[1]
        filename = os.path.join(directory, method, sample, f"rank_sum_weighted_out_{step_size}.tsv")
        if metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].isnull().values.any():
            continue
        if not os.path.isfile(filename):
            output_df.append([tumor, normal, 0, 0, 0, 0,
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])
            continue
        weighted_rank_sum_out = pd.read_csv(filename, sep="\t", header=0)
        if weighted_rank_sum_out.shape[0] == 0:
            output_df.append([tumor, normal, 0, 0, 0, 0,
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])
            continue
        features = weighted_rank_sum_out.loc[:, get_relevant_features_neofox()]
        for i, row in features.iterrows():
            if np.isnan(row['Selfsimilarity_conserved_binder']):
                features.loc[i, 'Selfsimilarity_conserved_binder'] = 0 if (features.loc[i, 'DAI'] > 0) else 1
        features = np.array(features)
        pred_proba = rf.xgb_model.predict_proba(features)[:, 1]
        output_df.append([tumor, normal, 
                          sum_top_x(weighted_rank_sum_out, 'qscore', topx), 
                          weighted_rank_sum_out[weighted_rank_sum_out['qscore'] > 0.54].shape[0],
                          sum_top_x(pd.DataFrame(pred_proba, columns=['pred']), 'pred', topx),
                          sum(1 for pred in pred_proba if pred > 0.6),
                          metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])


output_df = pd.DataFrame(output_df, columns=['TUMOR', 'NORMAL', f'QSCORE_MEDIAN_{topx}', 'LOAD_QSCORE', f'XGB_PRED_MEDIAN_{topx}', 'LOAD_XGB', 'RESPONSE'])
output_df["RESPONSE_BINARY"] = output_df["RESPONSE"].apply(lambda x: response_map[x])

recist_score = [response_map[r] for r in output_df['RESPONSE']]

# QScore as predictor for RECIST
pred_labels = (output_df[f'QSCORE_MEDIAN_{topx}'] > 0.54).astype(int)

accuracy = accuracy_score(recist_score, pred_labels)
precision = precision_score(recist_score, pred_labels)
recall = recall_score(recist_score, pred_labels)
f1 = f1_score(recist_score, pred_labels)
auc = roc_auc_score(recist_score, output_df[f'QSCORE_MEDIAN_{topx}'])

print("QScore Accuracy:", accuracy)
print("QScore Precision:", precision)
print("QScore Recall:", recall)
print("QScore F1-score:", f1)
print("QScore AUC:", auc)

stat, p_value = wilcoxon(output_df[f'QSCORE_MEDIAN_{topx}'], output_df['RESPONSE_BINARY'])
print(f"QScore Wilcoxon Statistic: {stat:.4f}, P-Value: {p_value:.4f}")

plt.figure(figsize=(4,6))
gfg = sns.boxplot(data=output_df, x="RESPONSE_BINARY", y=f'QSCORE_MEDIAN_{topx}', color='#94B6D2')
gfg.set(xlabel='RECIST classification', ylabel=f'top {topx} median immunogenicity score')
labels = gfg.get_xticklabels()
plt.ylim([0, 1])
gfg.set_xticklabels(['NR' if l.get_text() == '0' else 'R' for l in labels])
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/quality_score_analysis/images/responder_binary_top{topx}_qscore_{step_size}_boxplot.png")
plt.clf()

# Antigen load with QScore threshold as predictor for RECIST
best_f1 = 0
best_threshold = 0
for tr in range(1, 50):
    pred_labels = (output_df['LOAD_QSCORE'] > tr).astype(int)
    f1 = f1_score(recist_score, pred_labels)
    if f1 > best_f1:
        best_f1 = f1
        best_threshold = tr
print("Load QScore Best Threshold", best_threshold)

pred_labels = (output_df['LOAD_QSCORE'] > best_threshold).astype(int)

accuracy = accuracy_score(recist_score, pred_labels)
precision = precision_score(recist_score, pred_labels)
recall = recall_score(recist_score, pred_labels)
f1 = f1_score(recist_score, pred_labels)
auc = roc_auc_score(recist_score, output_df['LOAD_QSCORE'])

print("Load QScore Accuracy:", accuracy)
print("Load QScore Precision:", precision)
print("Load QScore Recall:", recall)
print("Load QScore F1-score:", f1)
print("Load QScore AUC:", auc)

stat, p_value = wilcoxon(output_df['LOAD_QSCORE'], output_df['RESPONSE_BINARY'])
print(f"Load QScore Wilcoxon Statistic: {stat:.4f}, P-Value: {p_value:.4f}")

plt.figure(figsize=(4,6))
gfg = sns.boxplot(data=output_df, x="RESPONSE_BINARY", y="LOAD_QSCORE", color='#94B6D2')
plt.ylim([0, 100])
gfg.set(xlabel='RECIST classification', ylabel=f'neoepitope load')
labels = gfg.get_xticklabels()
gfg.set_xticklabels(['NR' if l.get_text() == '0' else 'R' for l in labels])
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/quality_score_analysis/images/load_qscore_responder_binary_{step_size}_boxplot.png")
plt.clf()

plt.figure(figsize=(6,4))
sns.histplot(output_df[output_df['RESPONSE_BINARY'] == 0]['LOAD_QSCORE'], kde=True, bins=50, color="#94B6D2", binwidth=2, label="NR")
sns.histplot(output_df[output_df['RESPONSE_BINARY'] == 1]['LOAD_QSCORE'], kde=True, bins=50, color="#A5AB81", binwidth=2, label="R")
plt.xlim([0,100])
plt.legend()
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/quality_score_analysis/images/load_qscore_responder_binary_{step_size}_hist.png")
plt.clf()

# XGBoost prediction as predictor for RECIST
pred_labels = (output_df[f'XGB_PRED_MEDIAN_{topx}'] > 0.6).astype(int)

accuracy = accuracy_score(recist_score, pred_labels)
precision = precision_score(recist_score, pred_labels)
recall = recall_score(recist_score, pred_labels)
f1 = f1_score(recist_score, pred_labels)
auc = roc_auc_score(recist_score, output_df[f'XGB_PRED_MEDIAN_{topx}'])

print("XGB Accuracy:", accuracy)
print("XGB Precision:", precision)
print("XGB Recall:", recall)
print("XGB F1-score:", f1)
print("XGB AUC:", auc)

stat, p_value = wilcoxon(output_df[f'XGB_PRED_MEDIAN_{topx}'], output_df['RESPONSE_BINARY'])
print(f"XGB Wilcoxon Statistic: {stat:.4f}, P-Value: {p_value:.4f}")

plt.figure(figsize=(4,6))
gfg = sns.boxplot(data=output_df, x="RESPONSE_BINARY", y=f'XGB_PRED_MEDIAN_{topx}', color='#94B6D2')
gfg.set(xlabel='RECIST classification', ylabel=f'top {topx} median XGBoost pred. prob.')
labels = gfg.get_xticklabels()
gfg.set_xticklabels(['NR' if l.get_text() == '0' else 'R' for l in labels])
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/quality_score_analysis/images/xgb_responder_binary_top{topx}_qscore_{step_size}_boxplot")
plt.close()

# Antigen load with XGBoost prediction prob threshold as predictor for RECIST
best_f1 = 0
best_threshold = 0
for tr in range(1, 50):
    pred_labels = (output_df['LOAD_XGB'] > tr).astype(int)
    f1 = f1_score(recist_score, pred_labels)
    if f1 > best_f1:
        best_f1 = f1
        best_threshold = tr
print("Load XGB Best Threshold", best_threshold)
best_threshold = 10
pred_labels = (output_df['LOAD_XGB'] > best_threshold).astype(int)

accuracy = accuracy_score(recist_score, pred_labels)
precision = precision_score(recist_score, pred_labels)
recall = recall_score(recist_score, pred_labels)
f1 = f1_score(recist_score, pred_labels)
auc = roc_auc_score(recist_score, output_df['LOAD_XGB'])

print("Load XGB Accuracy:", accuracy)
print("Load XGB Precision:", precision)
print("Load XGB Recall:", recall)
print("Load XGB F1-score:", f1)
print("Load XGB AUC:", auc)

stat, p_value = wilcoxon(output_df['LOAD_XGB'], output_df['RESPONSE_BINARY'])
print(f"Load XGB Wilcoxon Statistic: {stat:.4f}, P-Value: {p_value:.4f}")

plt.figure(figsize=(4,6))
gfg = sns.boxplot(data=output_df, x="RESPONSE_BINARY", y="LOAD_XGB", color='#94B6D2')
plt.ylim([0, 100])
gfg.set(xlabel='RECIST classification', ylabel=f'neoepitope load')
labels = gfg.get_xticklabels()
gfg.set_xticklabels(['NR' if l.get_text() == '0' else 'R' for l in labels])
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/quality_score_analysis/images/load_xgb_responder_binary_{step_size}_boxplot.png")
plt.clf()

plt.figure(figsize=(6,4))
sns.histplot(output_df[output_df['RESPONSE_BINARY'] == 0]['LOAD_XGB'], kde=True, bins=50, color="#94B6D2", binwidth=2, label="NR")
sns.histplot(output_df[output_df['RESPONSE_BINARY'] == 1]['LOAD_XGB'], kde=True, bins=50, color="#A5AB81", binwidth=2, label="R")
plt.xlim([0,100])
plt.legend()
plt.tight_layout()
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/quality_score_analysis/images/load_xgb_responder_binary_{step_size}_hist.png")
plt.clf()