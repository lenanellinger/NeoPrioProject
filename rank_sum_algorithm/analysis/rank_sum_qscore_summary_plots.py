import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns

step_size = 25

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"
features_qscore = None
patient_id = 1
for cohort in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, cohort)):
        continue
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            if sample.startswith("_no"):
                continue
            tumor = sample.split("-")[0]
            normal = sample.split("-")[1]
            filename = os.path.join(directory, cohort, method, sample, f"rank_sum_weighted_out_{step_size}.tsv")
            if not os.path.isfile(filename):
                continue
            weighted_rank_sum_out = pd.read_csv(filename, sep="\t", header=0)
            weighted_rank_sum_out['patient_id'] = patient_id
            patient_id += 1
            if features_qscore is None:
                features_qscore = weighted_rank_sum_out
            else:
                features_qscore = pd.concat([features_qscore, weighted_rank_sum_out])
features_qscore.reset_index(drop=True, inplace=True)

# box plot
stats = features_qscore.groupby('patient_id').agg({'qscore': ['mean', 'min', 'max']})
stats.columns = ['mean', 'min', 'max']
stats.reset_index(inplace=True)

stats_melted = stats.melt(id_vars='patient_id', value_vars=['mean', 'min', 'max'],
                    var_name='Statistic', value_name='Value')

sns.boxplot(data=stats_melted, x="Statistic", y="Value")
plt.title('Distribution of Summary Statistics of QScore Across Patients')
plt.xlabel('Statistics')
plt.ylabel('Values')
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/qscore_summary_stats_{step_size}.png")