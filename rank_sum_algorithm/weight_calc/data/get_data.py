import os
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

sys.path.append(sys.path[0] + '/../../../analysis')
from helpers.get_data import get_relevant_features_neofox

directory_ml = "/mnt/storage2/users/ahnelll1/master_thesis/output_training_data"

def get_feature_data_NEPdb():
    """
    returns a DataFrame with neofox features
    """

    features = pd.read_csv(os.path.join(directory_ml,  "NEPdb_neofox_annotations.tsv"), sep="\t", header=0)
            
    return features


def get_train_data():
    features_with_meta = get_feature_data_NEPdb()

    labels = np.array(features_with_meta['response'])

    all_feature_names = get_relevant_features_neofox()

    features_with_meta = features_with_meta.loc[:, ['patientIdentifier'] + [name for name in all_feature_names if
                                                                            not name.startswith("Priority_score")]]
    features_with_meta['Selfsimilarity_conserved_binder'] = features_with_meta[
        'Selfsimilarity_conserved_binder'].fillna(0)

    features_with_meta_array = np.array(features_with_meta)

    train_features, test_features, train_labels, test_labels = train_test_split(features_with_meta_array, labels,
                                                                                test_size=0.2, random_state=42)

    train_features = train_features[:, 1:]
    train_labels_int = [0 if label == 'N' else 1 for label in train_labels]
    test_patients = test_features[:, 0]
    test_features = test_features[:, 1:]
    test_labels_int = [0 if label == 'N' else 1 for label in test_labels]

    imputer = SimpleImputer(strategy='median')
    imputer.fit(train_features)
    train_features = imputer.transform(train_features)
    test_features = imputer.transform(test_features)

    # TODO class imbalance

    return train_features, train_labels, train_labels_int, test_features, test_labels, test_labels_int, test_patients