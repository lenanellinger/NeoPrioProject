import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_squared_error

from scipy.stats import pearsonr

sys.path.append(sys.path[0] + '/../..')
from analysis.helpers.get_data import get_relevant_features_neofox
from rank_sum_algorithm.weight_calc.random_forest_model import RandomForestModel

rf = RandomForestModel()
step_size = 25

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"
features_qscore = None
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
            if features_qscore is None:
                features_qscore = weighted_rank_sum_out
            else:
                features_qscore = pd.concat([features_qscore, weighted_rank_sum_out])
features_qscore.reset_index(drop=True, inplace=True)
qscores = features_qscore.loc[:, ['qscore']]
features = features_qscore.loc[:, [name for name in get_relevant_features_neofox() if not name.startswith("Priority_score")]]

for i, row in features.iterrows():
    if np.isnan(row['Selfsimilarity_conserved_binder']):
        features.loc[i, 'Selfsimilarity_conserved_binder'] = 0 if (features.loc[i, 'DAI'] > 0) else 1
print(features['Selfsimilarity_conserved_binder'])
features = np.array(features)
qscores = np.array(qscores).flatten()

imputer = SimpleImputer(strategy='median')
imputer.fit(features)
features = imputer.transform(features)

pred_proba = rf.get_classifier().predict_proba(features)[:, 1]

# Correlation
pearson_corr = pearsonr(qscores, pred_proba)
print("Pearson Correlation", pearson_corr)

# Mean Square Error
mse = mean_squared_error(qscores, pred_proba)
print("Mean Square Error", mse)

# Scatter plot
plt.figure(figsize=(6, 6))
plt.scatter(qscores, pred_proba, alpha=0.5)
plt.plot([min(qscores), max(qscores)], [min(qscores), max(qscores)], color='red', linestyle='--')  # y=x line
plt.xlabel("Predictions 1")
plt.ylabel("Predictions 2")
plt.title(f"Scatter Plot (Pearson r = {pearson_corr.statistic:.2f})")
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/qscore_rf_scatter_{step_size}.png")

# Bland-Altman Plot
mean_pred = (qscores + pred_proba) / 2
diff_pred = qscores - pred_proba
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
plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/analysis/images/qscore_rf_bland_altman_{step_size}.png")