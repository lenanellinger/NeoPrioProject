import os
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import seaborn as sns
import matplotlib.pyplot as plt

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output/AxelMelanomaPhD"

metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
response_map = {
    'PD': 0,
    'SD': 0,
    'MR': 0,
    'CR': 1,
    'PR': 1
}

neoantigen_load = []
step_size = 25

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
            neoantigen_load.append([tumor, normal, 0, 
                metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])
            continue
        neoepitopes = pd.read_csv(filename, sep="\t", header=0)
        neoantigen_load.append([tumor, normal, neoepitopes.shape[0], metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]])

df_neoantigen_load = pd.DataFrame(neoantigen_load, columns=['TUMOR', 'NORMAL', 'LOAD', 'RESPONSE'])
df_neoantigen_load["RESPONSE_BINARY"] = df_neoantigen_load["RESPONSE"].apply(lambda x: response_map[x])

recist_score = [response_map[r] for r in df_neoantigen_load['RESPONSE']]

# Boxplot
sns.boxplot(data=df_neoantigen_load, x="RESPONSE", y="LOAD")
plt.ylim([0, 500])
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/load_responder.png")
plt.clf()
sns.boxplot(data=df_neoantigen_load, x="RESPONSE_BINARY", y="LOAD")
plt.ylim([0, 500])
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/load_responder_binary.png")