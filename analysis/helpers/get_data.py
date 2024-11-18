import os
import pandas as pd

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

def get_feature_data(prefilter=True, cohorts=['AxelMelanomaPhD', 'SomaticAndTreatment']):
    """
    returns a DataFrame (filtered by binding affinity prefilter) for given cohorts
    """
    features_data = None

    for cohort in cohorts:
        file = os.path.join(directory, cohort, cohort + "_all.tsv")
        features = pd.read_csv(file, sep="\t", header=0)
        
        weights = {
            'MHCflurry MT IC50 Score': 0.4,
            'NetMHCpan MT IC50 Score': 0.4,
            'NetMHC MT IC50 Score': 0.2
        }
        
        features['IC50 mean'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].mean(axis=1)
        
        def weighted_mean(row):
            score_sum = 0
            weight_sum = 0
            for feature, weight in weights.items():
                if not pd.isna(row[feature]):  # Skip NaN values
                    score_sum += row[feature] * weight
                    weight_sum += weight
            return score_sum / weight_sum if weight_sum > 0 else np.nan

        features['IC50 weighted mean'] = features.apply(weighted_mean, axis=1)
        features['cohort'] = cohort
        
        if prefilter:
            features = features[features['IC50 weighted mean'] < 500]

        if features_data is None:
            features_data = features
        else:
            features_data = pd.concat([features_data, features], axis=0, ignore_index=True)
            
    return features_data
    
def get_relevant_features(details=False):
    """
    returns relevant feature column names (or as detailed object)
    """
    return get_relevant_features_neofox(details) +  get_relevant_features_pvacseq(details)
    
    
def get_relevant_features_neofox(details=False):
    """
    returns relevant neofox feature column names (or as detailed object)
    """
    features = [
        {
        #    'name': 'amplitude', 
        #    'quantile': 'upper',
        #    'interval': [0, 10],
        #    'cutoff': {
        #        'use_percentage': False,
        #        'good_distribution': True,
        #        'specific_value': 1
        #    }
        #},{
            'name': 'DAI', 
            'quantile': 'upper',
            'interval': [-915, 3000],
            'cutoff': {
                'use_percentage': False,
                'good_distribution': True,
                'specific_value': 0
            }
        },{
            'name': 'dissimilarity_score', 
            'quantile': 'upper'
        },{
            'name': 'IEDB_Immunogenicity', 
            'quantile': 'upper',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': True
            }
        },{
            'name': 'MixMHCpred_score', 
            'quantile': 'upper',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': True
            }
        #},{
        #    'name': 'pathogen_similarity', 
        #    'quantile': 'upper',
        #    'interval': [0.99, 1]
        },{
            'name': 'Selfsimilarity_conserved_binder', 
            'quantile': 'lower',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': False
            }
        },{
            'name': 'recognition_potential', 
            'quantile': 'upper'
        },{
            'name': 'Priority_score_imputed_fromDNA', 
            'quantile': 'upper',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': False
            }
        },{
            'name': 'Priority_score_imputed_fromRNA', 
            'quantile': 'upper',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': False
            }
        #},{
        #    'name': 'Tcell_predictor', 
        #    'quantile': 'upper',
        #    'cutoff': {
        #        'use_percentage': True,
        #        'good_distribution': False
        #    }
        },{	
            'name': 'hex_alignment_score', 
            'quantile': 'upper',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': True
            }
        },{
            'name': 'PRIME_score', 
            'quantile': 'upper',
            'cutoff': {
                'use_percentage': True,
                'good_distribution': True
            }
        #},{
        #    'name': 'Predicted Stability',
        #    'quantile': 'upper',
        #    'cutoff': {
        #        'use_percentage': True,
        #        'good_distribution': True
        #    }
        },{
            'name': 'imputedGeneExpression', 
            'quantile': 'upper',
            'interval': [0, 40],
            'cutoff': {
                'use_percentage': True,
                'good_distribution': True
            }
        }
    ]
    
    if details:
        return features
    else:
        return [f['name'] for f in features]
    
def get_relevant_features_pvacseq(details=False):
    """
    returns relevant pvacseq feature column names (or as detailed object)
    """
    features = [
        {
            'name': 'NetMHC MT IC50 Score', 
            'quantile': 'lower'
        },{
            'name': 'NetMHCpan MT IC50 Score', 
            'quantile': 'lower'
        },{
            'name': 'MHCflurry MT IC50 Score', 
            'quantile': 'lower'
        }
    ]
    
    if details:
        return features
    else:
        return [f['name'] for f in features]