import pandas as pd
import numpy as np

rf_output = pd.read_csv('/mnt/storage2/users/ahnelll1/master_thesis/training_data/NEPdb_neofox_annotations_random_forest_out.tsv', sep='\t', header=0)
rs_output = pd.read_csv('/mnt/storage2/users/ahnelll1/master_thesis/training_data/NEPdb_neofox_annotations_rank_sum_out_25.tsv', sep='\t', header=0)

comparison_df = rf_output.set_index('patientIdentifier').join(rs_output.set_index('patientIdentifier')).loc[:, ['qscore', 'pred_proba']]
comparison_df['diff_score'] = np.abs(comparison_df['qscore'] - comparison_df['pred_proba'])


threshold = 0.8
comparison_df['qscore_class'] = comparison_df['qscore'] > threshold
comparison_df['pred_proba_class'] = comparison_df['pred_proba'] > threshold
comparison_df['diff_class'] = (comparison_df['qscore_class'] & ~comparison_df['pred_proba_class']) | (~comparison_df['qscore_class'] & comparison_df['pred_proba_class'])

print(comparison_df)
print('Deviation Score:', comparison_df.mean()['diff_score'])
print('Deviation Class:', comparison_df['diff_class'].sum() / comparison_df.shape[0])