import os
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import seaborn as sns
import matplotlib.pyplot as plt

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

def sum_top_10(df, column):
    sorted_scores = df[column].sort_values(ascending=False)
    top_10_scores = sorted_scores.head(10)
    
    return top_10_scores.sum()

for method in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, method)):
        continue
    for sample in os.listdir(os.path.join(directory, method)):
        if sample.startswith("_no"):
            continue
        tumor = sample.split("-")[0]
        normal = sample.split("-")[1]
        filename = os.path.join(directory, method, sample, "rank_sum_weighted_out_25.tsv")
        if metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].isnull().values.any():
            continue
        if not os.path.isfile(filename):
            prediction_weighted.append([tumor, normal, 0, 
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])
            continue
        weighted_rank_sum_out = pd.read_csv(filename, sep="\t", header=0)
        if weighted_rank_sum_out.shape[0] == 0:
            prediction_weighted.append([tumor, normal, 0, 
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])
            continue
        prediction_weighted.append([tumor, normal, sum_top_10(weighted_rank_sum_out, 'qscore'), metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])

df_prediction_weighted = pd.DataFrame(prediction_weighted, columns=['TUMOR', 'NORMAL', 'QSCORE', 'RESPONSE'])
pred_labels = (df_prediction_weighted['QSCORE'] > 5.5).astype(int)

recist_score = [response_map[r] for r in df_prediction_weighted['RESPONSE']]

accuracy = accuracy_score(recist_score, pred_labels)
precision = precision_score(recist_score, pred_labels)
recall = recall_score(recist_score, pred_labels)
f1 = f1_score(recist_score, pred_labels)
auc = roc_auc_score(recist_score, df_prediction_weighted['QSCORE'])

print("Accuracy:", accuracy)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc)

# Boxplot
sns.boxplot(data=df_prediction_weighted, x="RESPONSE", y="QSCORE")
plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/images/responder_top10_qscore_boxplot_v1.png")  