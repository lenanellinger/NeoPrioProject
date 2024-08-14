import os
import pandas as pd
import numpy as np
import sys
from optparse import OptionParser
import json
import bisect


def calc_qscore(row, columns, number_of_quantiles):
    qscore_series = row[[s for s in columns if s.endswith('_qscore')]]
    return qscore_series.sum(skipna=True) / (number_of_quantiles * (~qscore_series.isnull()).sum())
    
def rank_sum(pvacseq_filename, neofox_filename):  
    pvacseq_df = pd.read_csv(pvacseq_filename, sep="\t", header=0)
    neofox_df = pd.read_csv(neofox_filename, sep="\t", header=0)
    
    f = open('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/data/quantiles_25.json')
    features = json.load(f)
    f.close()
    
    if neofox_df.shape[0] != pvacseq_df.shape[0]:
        raise ValueError("Neofox output and pVACseq output should have the same amount of rows.")
    annotation_df = neofox_df.merge(right=pvacseq_df.loc[:, ['Chromosome', 'Start', 'Stop', 'Transcript', 'Variant Type', 
    'MHCflurry WT IC50 Score', 'MHCflurry MT IC50 Score', 'MHCflurry WT Percentile', 'MHCflurry MT Percentile', 
    'NetMHC WT IC50 Score', 'NetMHC MT IC50 Score', 'NetMHC WT Percentile', 'NetMHC MT Percentile', 
    'NetMHCpan WT IC50 Score', 'NetMHCpan MT IC50 Score', 'NetMHCpan WT Percentile',  'NetMHCpan MT Percentile', 
    'Predicted Stability', 'Half Life', 'Stability Rank', 'NetMHCstab allele']], on=['Chromosome', 'Start', 'Stop', 'Transcript'])
    
    # Prefilter (majority vote on binding affinities below 500)
    voting = pd.DataFrame()
    voting['NetMHC'] = annotation_df['NetMHC MT IC50 Score'] < 500
    voting['NetMHCpan'] = annotation_df['NetMHCpan MT IC50 Score'] < 500
    voting['MHCflurry'] = annotation_df['MHCflurry MT IC50 Score'] < 500
    voting['majority'] = voting.mode(axis=1)[0]
    annotation_df_prefiltered = annotation_df[voting['majority']]
    
    output_df = annotation_df_prefiltered.loc[:, ['Transcript', 'Chromosome', 'Start', 'Stop'] + [f for f in features]]
    number_of_quantiles = None
    
    for feature_name in features:
        index = output_df.columns.get_loc(feature_name)
        number_of_quantiles = len(features[feature_name]['quantiles'])
        
        # Ranking
        column_name_rank = feature_name + '_rank'
        if feature_name == 'Selfsimilarity_conserved_binder':
            # if Selfsimilarity is NaN, the DAI is almost zero -> NaN would be preferable
            rank_column = output_df[feature_name].rank(method='max', ascending=features[feature_name]['direction'] == 'lower', na_option='top')
        elif feature_name == 'imputedGeneExpression':
            output_df.insert(index, 'rnaExpression', annotation_df_prefiltered['rnaExpression'])
            index += 1
            # if RNA data given then measured RNA expression should be used, else the imputed
            rank_column = output_df['rnaExpression'].where(~output_df['rnaExpression'].isnull(), output_df['imputedGeneExpression']).rank(method='max', ascending=features[feature_name]['direction'] == 'lower', na_option='bottom')             
            column_name_rank = 'gene_expression_rank'
        elif feature_name == 'Priority_score_imputed_fromDNA' or feature_name == 'Priority_score_imputed_fromRNA':
            output_df.insert(index, feature_name.replace('_imputed', ''), annotation_df_prefiltered[feature_name.replace('_imputed', '')])
            index += 1
            # if RNA data given then measured RNA expression should be used for priority score, else the imputed            
            rank_column = output_df[feature_name.replace('_imputed', '')].where(~output_df[feature_name.replace('_imputed', '')].isnull(), output_df[feature_name]).rank(method='max', ascending=features[feature_name]['direction'] == 'lower', na_option='bottom')            
            column_name_rank = feature_name.replace('_imputed', '') + '_rank'
        else:
            rank_column = output_df[feature_name].rank(method='max', ascending=features[feature_name]['direction'] == 'lower', na_option='bottom')
        output_df.insert(index + 1, column_name_rank, rank_column)
            
        # Quantile Score
        column_name_qscore = feature_name + '_qscore'
        if feature_name == 'imputedGeneExpression':
            # if RNA data given then measured RNA expression should be used, else the imputed
            quantile_score_column = [bisect.bisect_left(features[feature_name]['quantiles'], num) if not np.isnan(num) else None for _, num in output_df['rnaExpression'].where(~output_df['rnaExpression'].isnull(), output_df['imputedGeneExpression']).items()]
            column_name_qscore = 'gene_expression_qscore'
        elif feature_name == 'Priority_score_imputed_fromDNA' or feature_name == 'Priority_score_imputed_fromRNA':
            # if RNA data given then measured RNA expression should be used for priority score, else the imputed  
            quantile_score_column = [bisect.bisect_left(features[feature_name]['quantiles'], num) if not np.isnan(num) else None for _, num in output_df[feature_name.replace('_imputed', '')].where(~output_df[feature_name.replace('_imputed', '')].isnull(), output_df[feature_name]).items()]            
            column_name_qscore = feature_name.replace('_imputed', '') + '_qscore'
        else:
            quantile_score_column = [bisect.bisect_left(features[feature_name]['quantiles'], num) if not np.isnan(num) else None for _, num in output_df[feature_name].items()]
        quantile_score_column = quantile_score_column if features[feature_name]['direction'] == 'upper' else [number_of_quantiles - x if not pd.isnull(x) else None for x in quantile_score_column]
        output_df.insert(index + 2, column_name_qscore, quantile_score_column)
      
    # Output  
    output_df['rank_sum'] = output_df[[s for s in list(output_df.columns) if s.endswith('_rank')]].sum(axis=1, skipna=False)
    output_df['qscore'] = None
    if not output_df.shape[0] == 0:
        output_df['qscore'] = output_df.apply(calc_qscore, axis='columns', args=(list(output_df.columns), number_of_quantiles))

    output_df = output_df.sort_values(by='rank_sum')
            
    output_df.to_csv(os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(neofox_filename))), 'rank_sum_out.tsv'), sep='\t', index=False, header=True)

def main(argv):
    usage = "usage: python neoantigen_prioritization.py "
    desc = "Prioritizes neoantigens based on neofox features with a rank sum algorithm."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--pvacseq-filtered-output", action="store", dest="pvacseq_filtered_output", help="Filtered TSV file output of pVACseq")
    parser.add_option("--neofox-output", action="store", dest="neofox_output", help="TSV file output of neofox")
    (options, args) = parser.parse_args()
    
    rank_sum(options.pvacseq_filtered_output, options.neofox_output)

    
if __name__ == "__main__":
    main(sys.argv)