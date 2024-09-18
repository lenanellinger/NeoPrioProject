import os
import pandas as pd
import numpy as np
import sys
from optparse import OptionParser
import json
import bisect


def calc_qscore(row, columns, number_of_quantiles, weights):
    qscore_series = row[[s for s in columns if s.endswith('_qscore')]]
    qscore_series_renamed = qscore_series.rename(index=lambda x: x.replace('_qscore', ''))
    valid_scores = qscore_series_renamed.dropna()

    valid_weights = {k: weights[k] for k in valid_scores.index}
    total_valid_weight = sum(valid_weights.values())
    adjusted_weights = {k: v / total_valid_weight for k, v in valid_weights.items()}

    weighted_qscores = (valid_scores / number_of_quantiles).mul(pd.Series(adjusted_weights))
    return weighted_qscores.sum(skipna=True)


def create_annotations_df(pvacseq_filename, neofox_filename):
    pvacseq_df = pd.read_csv(pvacseq_filename, sep="\t", header=0)
    neofox_df = pd.read_csv(neofox_filename, sep="\t", header=0)

    if neofox_df.shape[0] != pvacseq_df.shape[0]:
        raise ValueError("Neofox output and pVACseq output should have the same amount of rows.")
    annotation_df = neofox_df.merge(
        right=pvacseq_df.loc[:, ['Chromosome', 'Start', 'Stop', 'Transcript', 'Variant Type',
                                 'MHCflurry WT IC50 Score', 'MHCflurry MT IC50 Score', 'MHCflurry WT Percentile',
                                 'MHCflurry MT Percentile',
                                 'NetMHC WT IC50 Score', 'NetMHC MT IC50 Score', 'NetMHC WT Percentile',
                                 'NetMHC MT Percentile',
                                 'NetMHCpan WT IC50 Score', 'NetMHCpan MT IC50 Score', 'NetMHCpan WT Percentile',
                                 'NetMHCpan MT Percentile',
                                 'Predicted Stability', 'Half Life', 'Stability Rank', 'NetMHCstab allele']],
        on=['Chromosome', 'Start', 'Stop', 'Transcript'])

    # Prefilter (weighted average below 500)
    annotation_df_prefiltered = annotation_df[(annotation_df['MHCflurry MT IC50 Score'] * 0.4 + annotation_df[
        'NetMHCpan MT IC50 Score'] * 0.4 + annotation_df['NetMHC MT IC50 Score'] * 0.2) < 500]

    return annotation_df_prefiltered


def rank_sum_qscore(pvacseq_filename, neofox_filename, step_size):
    output_df = ranking(create_annotations_df(pvacseq_filename, neofox_filename),
                        ['Transcript', 'Chromosome', 'Start', 'Stop'], step_size)

    output_df.to_csv(
        os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(neofox_filename))),
                     f'rank_sum_weighted_out_{step_size}.tsv'),
        sep='\t', index=False, header=True)


def rank_sum_qscore_training_data(neofox_filename, step_size):
    neofox_df = pd.read_csv(neofox_filename, sep="\t", header=0)

    output_df = ranking(neofox_df, ['patientIdentifier'], step_size)

    output_df.to_csv(os.path.join(os.path.abspath(os.path.dirname(neofox_filename)),
                                  os.path.splitext(os.path.basename(neofox_filename))[
                                      0] + f'_rank_sum_weighted_out_{step_size}.tsv'),
                     sep='\t', index=False, header=True)


def ranking(input_df, meta_columns, step_size):
    f = open(f'/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/data/quantiles_{step_size}.json')
    features = json.load(f)
    f.close()

    f = open('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/weight_calc/weights.json')
    weights = json.load(f)
    f.close()

    output_df = input_df.loc[:, meta_columns + [f for f in features if f in input_df]]

    number_of_quantiles = None
    for feature_name in features:
        if feature_name not in output_df.columns:
            continue
        index = output_df.columns.get_loc(feature_name)
        number_of_quantiles = len(features[feature_name]['quantiles'])

        # Ranking
        column_name_rank = feature_name + '_rank'
        if feature_name == 'Selfsimilarity_conserved_binder':
            # if Selfsimilarity is NaN, the DAI is almost zero -> NaN would be preferable
            rank_column = output_df[feature_name].rank(method='max',
                                                       ascending=features[feature_name]['direction'] == 'lower',
                                                       na_option='top')
        elif feature_name == 'imputedGeneExpression':
            output_df.insert(index, 'rnaExpression', input_df['rnaExpression'])
            index += 1
            # if RNA data given then measured RNA expression should be used, else the imputed
            rank_column = output_df['rnaExpression'].where(~output_df['rnaExpression'].isnull(),
                                                           output_df['imputedGeneExpression']).rank(method='max',
                                                                                                    ascending=features[
                                                                                                                  feature_name][
                                                                                                                  'direction'] == 'lower',
                                                                                                    na_option='bottom')
            column_name_rank = 'gene_expression_rank'
        elif feature_name == 'Priority_score_imputed_fromDNA' or feature_name == 'Priority_score_imputed_fromRNA':
            output_df.insert(index, feature_name.replace('_imputed', ''),
                             input_df[feature_name.replace('_imputed', '')])
            index += 1
            # if RNA data given then measured RNA expression should be used for priority score, else the imputed            
            rank_column = output_df[feature_name.replace('_imputed', '')].where(
                ~output_df[feature_name.replace('_imputed', '')].isnull(), output_df[feature_name]).rank(method='max',
                                                                                                         ascending=
                                                                                                         features[
                                                                                                             feature_name][
                                                                                                             'direction'] == 'lower',
                                                                                                         na_option='bottom')
            column_name_rank = feature_name.replace('_imputed', '') + '_rank'
        else:
            rank_column = output_df[feature_name].rank(method='max',
                                                       ascending=features[feature_name]['direction'] == 'lower',
                                                       na_option='bottom')
        output_df.insert(index + 1, column_name_rank, rank_column)

        # Quantile Score
        column_name_qscore = feature_name + '_qscore'
        if feature_name == 'imputedGeneExpression':
            # if RNA data given then measured RNA expression should be used, else the imputed
            quantile_score_column = [
                bisect.bisect_left(features[feature_name]['quantiles'], num) if not np.isnan(num) else None for _, num
                in output_df['rnaExpression'].where(~output_df['rnaExpression'].isnull(),
                                                    output_df['imputedGeneExpression']).items()]
            column_name_qscore = 'gene_expression_qscore'
        elif feature_name == 'Priority_score_imputed_fromDNA' or feature_name == 'Priority_score_imputed_fromRNA':
            # if RNA data given then measured RNA expression should be used for priority score, else the imputed  
            quantile_score_column = [
                bisect.bisect_left(features[feature_name]['quantiles'], num) if not np.isnan(num) else None for _, num
                in output_df[feature_name.replace('_imputed', '')].where(
                    ~output_df[feature_name.replace('_imputed', '')].isnull(), output_df[feature_name]).items()]
            column_name_qscore = feature_name.replace('_imputed', '') + '_qscore'
        else:
            quantile_score_column = [
                bisect.bisect_left(features[feature_name]['quantiles'], num) if not np.isnan(num) else None for _, num
                in output_df[feature_name].items()]
        quantile_score_column = quantile_score_column if features[feature_name]['direction'] == 'upper' else [
            number_of_quantiles - x if not pd.isnull(x) else None for x in quantile_score_column]
        output_df.insert(index + 2, column_name_qscore, quantile_score_column)

    # Output
    ranks = output_df[[s for s in list(output_df.columns) if s.endswith('_rank')]]
    ranks_renamed = ranks.rename(columns=lambda x: x.replace('_rank', ''))

    valid_weights = {k: 1/weights[k] for k in weights} # use reciprocal of weights, because better rank is lower
    total_valid_weight = sum(valid_weights.values())
    adjusted_weights = {k: v / total_valid_weight for k, v in valid_weights.items()}

    weighted_ranks = ranks_renamed.mul(pd.Series(adjusted_weights),
                                       axis=1)
    output_df['rank_sum'] = weighted_ranks.sum(axis=1, skipna=False)
    output_df['qscore'] = None
    if not output_df.shape[0] == 0:
        output_df['qscore'] = output_df.apply(calc_qscore, axis='columns',
                                              args=(list(output_df.columns), number_of_quantiles, weights))

    output_df = output_df.sort_values(by='rank_sum')

    return output_df


def main(argv):
    usage = "usage: python neoantigen_prioritization.py "
    desc = "Prioritizes neoantigens based on neofox features with a rank sum algorithm."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--pvacseq-filtered-output", action="store", dest="pvacseq_filtered_output",
                      help="Filtered TSV file output of pVACseq")
    parser.add_option("--neofox-output", action="store", dest="neofox_output", help="TSV file output of neofox")
    (options, args) = parser.parse_args()

    rank_sum_qscore(options.pvacseq_filtered_output, options.neofox_output)


if __name__ == "__main__":
    main(sys.argv)
