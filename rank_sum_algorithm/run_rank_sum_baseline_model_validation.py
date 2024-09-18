import os
import pandas as pd
from neoantigen_prioritization_rank_sum import ranking, create_annotations_df
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

# create baseline
cohort = "SomaticAndTreatment"
baseline = None
for method in os.listdir(os.path.join(directory, cohort)):
    if os.path.isfile(os.path.join(directory, cohort, method)):
        continue
    for sample in os.listdir(os.path.join(directory, cohort, method)):
        if sample.startswith("_no"):
            continue
        pvacseq_filename = os.path.join(directory, cohort, method, sample, "pVACseq", "MHC_Class_I",
                                        sample.split("-")[0] + ".filtered.tsv")
        neofox_filename = os.path.join(directory, cohort, method, sample, "neofox",
                                       sample.split("-")[0] + "_neofox_annotations.tsv")
        if not os.path.isfile(neofox_filename):
            continue
        annotations_df = create_annotations_df(pvacseq_filename, neofox_filename)
        if baseline is None:
            baseline = annotations_df
        else:
            baseline = pd.concat([baseline, annotations_df])
baseline['cohort'] = 'baseline'

# validate with external cohort
metadata = pd.read_csv("/mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv", sep="\t", header=0)
response_map = {
    'PD': 0,
    'SD': 0,
    'MR': 0,
    'CR': 1,
    'PR': 1
}

cohort = "AxelMelanomaPhD"
predictions = []
responses = []
for method in os.listdir(os.path.join(directory, cohort)):
    if os.path.isfile(os.path.join(directory, cohort, method)):
        continue
    for sample in os.listdir(os.path.join(directory, cohort, method)):
        if sample.startswith("_no"):
            continue
        tumor = sample.split("-")[0]
        normal = sample.split("-")[1]
        if metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].isnull().values.any():
            continue
        pvacseq_filename = os.path.join(directory, cohort, method, sample, "pVACseq", "MHC_Class_I",
                                        sample.split("-")[0] + ".filtered.tsv")
        neofox_filename = os.path.join(directory, cohort, method, sample, "neofox",
                                       sample.split("-")[0] + "_neofox_annotations.tsv")
        if not os.path.isfile(neofox_filename):
            continue
        annotations_df = create_annotations_df(pvacseq_filename, neofox_filename)
        annotations_df['cohort'] = 'external'

        df_concat = pd.concat([baseline, annotations_df])

        output = ranking(df_concat, ['cohort'], 25)
        output_sorted = output.sort_values('rank_sum', ascending=True).reset_index(drop=True)
        response = response_map[metadata[(metadata['TUMOR'] == tumor) & (metadata['NORMAL'] == normal)]['RESPONSE'].values[0]]
        if output_sorted[output_sorted['cohort'] == 'external'].shape[0] > 0:
            position = output_sorted[output_sorted['cohort'] == 'external'].index[0]/output_sorted.shape[0]
        else:
            position = 1
        predictions.append(1-position)
        responses.append(response)

pred_labels = [1 if p > 0.99 else 0 for p in predictions]

accuracy = accuracy_score(responses, pred_labels)
precision = precision_score(responses, pred_labels)
recall = recall_score(responses, pred_labels)
f1 = f1_score(responses, pred_labels)
auc = roc_auc_score(responses, predictions)

print("Accuracy:", accuracy)
print("Precision:", precision)
print("Recall:", recall)
print("F1-score:", f1)
print("AUC:", auc)