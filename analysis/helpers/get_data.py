import os
import pandas as pd

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background"

def get_feature_data(majority_vote=True, cohorts=['AxelMelanomaPhD', 'SomaticAndTreatment']):
    """
    returns a DataFrame (filtered by majority vote of binding affinity) for given cohorts
    """
    features_data = None

    for cohort in cohorts:
        file = os.path.join(directory, cohort, cohort + "_all.tsv")
        features = pd.read_csv(file, sep="\t", header=0)
        
        #features['IC50 mean'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].mean(axis=1)
        #features['IC50 median'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].median(axis=1)
        features['cohort'] = cohort
        
        if majority_vote:
            voting = pd.DataFrame()
            voting['NetMHC'] = features['NetMHC MT IC50 Score'] < 500
            voting['NetMHCpan'] = features['NetMHCpan MT IC50 Score'] < 500
            voting['MHCflurry'] = features['MHCflurry MT IC50 Score'] < 500
            voting['majority'] = voting.mode(axis=1)[0]
            features = features[voting['majority']]

        if features_data is None:
            features_data = features
        else:
            features_data = pd.concat([features_data, features], axis=0, ignore_index=True)
            
    return features_data
    
def get_relevant_features(details=False):
    """
    returns relevant feature column names (as detailed obejct)
    """
    return get_relevant_features_neofox(details) +  get_relevant_features_pvacseq(details)
    
    
def get_relevant_features_neofox(details=False):
    """
    returns relevant neofox feature column names (as detailed obejct)
    """
    features = [
        {
            'name': 'amplitude', 
            'quantile': 'upper',
            'interval': [0, 10]
        },{
            'name': 'DAI', 
            'quantile': 'upper',
        },{
            'name': 'dissimilarity_score', 
            'quantile': 'upper',
            'interval': [0, 0.001]
        },{
            'name': 'IEDB_Immunogenicity', 
            'quantile': 'upper',
        },{
            'name': 'MixMHCpred_score', 
            'quantile': 'lower',
        },{
            'name': 'pathogen_similarity', 
            'quantile': 'upper',
            'interval': [0.99, 1]
        },{
            'name': 'Selfsimilarity_conserved_binder', 
            'quantile': 'lower',
        },{
            'name': 'recognition_potential', 
            'quantile': 'upper',
            'interval': [0, 0.54]
        },{
            'name': 'Priority_score_fromDNA', 
            'quantile': 'upper'
        },{
            'name': 'Priority_score_fromRNA', 
            'quantile': 'upper'
        },{
            'name': 'Priority_score_imputed_fromDNA', 
            'quantile': 'upper'
        },{
            'name': 'Priority_score_imputed_fromRNA', 
            'quantile': 'upper'
        },{
            'name': 'Tcell_predictor', 
            'quantile': 'upper'
        },{	
            'name': 'hex_alignment_score', 
            'quantile': 'upper'
        #},{
        #    'name': 'PRIME_score', 
        #    'quantile': 'upper'
        },{
            'name': 'Predicted Stability', 
            'quantile': 'upper'
        }
    ]
    
    if details:
        return features
    else:
        return [f['name'] for f in features]
    
def get_relevant_features_pvacseq(details=False):
    """
    returns relevant pvacseq feature column names (as detailed obejct)
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