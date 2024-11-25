import os
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

from imblearn.under_sampling import RandomUnderSampler, TomekLinks
from imblearn.over_sampling import RandomOverSampler, SMOTE

sys.path.append(sys.path[0] + '/../../../feature_analysis')
from helpers.get_data import get_relevant_features_neofox

directory_ml = "/mnt/storage2/users/ahnelll1/master_thesis/output_training_data"

def get_feature_data_NEPdb():
    """
    returns a DataFrame with neofox features
    """

    features = pd.read_csv(os.path.join(directory_ml,  "NEPdb_filtered_neofox.tsv"), sep="\t", header=0)
            
    return features

def plot_pca(name, x, y):
    pca = PCA(n_components=2)
    X = pca.fit_transform(x)
    colors = ['#1F77B4', '#FF7F0E']
    markers = ['o', 's']
    for l, c, m in zip(np.unique(y), colors, markers):
        plt.scatter(
            X[y==l, 0],
            X[y==l, 1],
            c=c, label=l, marker=m
        )
    plt.title('Imbalanced dataset (2 PCA components)')
    plt.legend(loc='upper right')
    plt.savefig(f"/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/data_pca_{name}.png")
    plt.clf()

def get_train_data():
    features_with_meta = get_feature_data_NEPdb()
    features_with_meta = features_with_meta[features_with_meta['affinityMutated'] <= 500]

    labels = np.array(features_with_meta['response'])

    all_feature_names = get_relevant_features_neofox()
    
    print("Number of Negative samples:", len([i for i in labels if i == 'N']))
    print("Number of Positive samples:", len([i for i in labels if i == 'P']))

    features_with_meta = features_with_meta.loc[:, ['patientIdentifier'] + all_feature_names]
    for i, row in features_with_meta.iterrows():
        if features_with_meta.isnull().loc[i, 'Selfsimilarity_conserved_binder']:
            features_with_meta.loc[i, 'Selfsimilarity_conserved_binder'] = 0 if (features_with_meta.loc[i, 'DAI'] > 0) else 1

    features_with_meta_array = np.array(features_with_meta)
    
    train_features, test_features_tmp, train_labels, test_labels_tmp = train_test_split(features_with_meta_array, labels,
                                                                                test_size=0.3, random_state=42, stratify=labels)

    val_features, test_features, val_labels, test_labels = train_test_split(test_features_tmp, test_labels_tmp,
                                                                                test_size=0.5, random_state=42, stratify=test_labels_tmp)
                                                                                
    train_features = train_features[:, 1:]
    train_labels_int = [0 if label == 'N' else 1 for label in train_labels]
    val_features = val_features[:, 1:]
    val_labels_int = [0 if label == 'N' else 1 for label in val_labels]
    test_patients = test_features[:, 0]
    test_features = test_features[:, 1:]
    test_labels_int = [0 if label == 'N' else 1 for label in test_labels]
    
    imputer = SimpleImputer(strategy='median')
    imputer.fit(train_features)
    train_features = imputer.transform(train_features)
    val_features = imputer.transform(val_features)
    test_features = imputer.transform(test_features)
    
    rus = SMOTE(sampling_strategy=0.5, random_state=42)
    train_features_res, train_labels_res = rus.fit_resample(train_features, train_labels)
    train_labels_int_res = [0 if label == 'N' else 1 for label in train_labels_res]
    
    print("After resample: Number of Negative samples:", len([i for i in train_labels_res if i == 'N']))
    print("After resample: Number of Positive samples:", len([i for i in train_labels_res if i == 'P']))
    
    plot_pca('dataset_train', train_features, train_labels)
    plot_pca('dataset_test', test_features, test_labels)

    return train_features, train_labels_int, train_features_res, train_labels_int_res, val_features, val_labels_int, test_features, test_labels_int, test_patients