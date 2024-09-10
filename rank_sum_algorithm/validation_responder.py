import os
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background/AxelMelanomaPhD"

metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
response_map = {
    'PD': 0,
    'SD': 0,
    'MR': 0,
    'CR': 1,
    'PR': 1
}

prediction_weighted = []

for method in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, method)):
        continue
    for sample in os.listdir(os.path.join(directory, method)):
        if sample.startswith("_no"):
            continue
        tumor = sample.split("-")[0]
        normal = sample.split("-")[1]
        filename = os.path.join(directory, method, sample, "rank_sum_weighted_out_25_v3.tsv")
        if metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].isnull().values.any():
            continue
        if not os.path.isfile(filename):
            prediction_weighted.append([tumor, normal, 0, response_map[
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]]])
            continue
        weighted_rank_sum_out = pd.read_csv(filename, sep="\t", header=0)
        if weighted_rank_sum_out.shape[0] == 0:
            prediction_weighted.append([tumor, normal, 0, response_map[
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]]])
            continue
        prediction_weighted.append([tumor, normal, weighted_rank_sum_out.max()['qscore'], response_map[metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]]])

df_prediction_weighted = pd.DataFrame(prediction_weighted, columns=['TUMOR', 'NORMAL', 'QSCORE', 'RESPONSE'])
print(df_prediction_weighted.to_string())
pred_labels = (df_prediction_weighted['QSCORE'] > 0.8).astype(int)

accuracy = accuracy_score(df_prediction_weighted['RESPONSE'], pred_labels)
precision = precision_score(df_prediction_weighted['RESPONSE'], pred_labels)
recall = recall_score(df_prediction_weighted['RESPONSE'], pred_labels)
f1 = f1_score(df_prediction_weighted['RESPONSE'], pred_labels)
auc = roc_auc_score(df_prediction_weighted['RESPONSE'], df_prediction_weighted['QSCORE'])

print("Accuracy:", accuracy)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc)