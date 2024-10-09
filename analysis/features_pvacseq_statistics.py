import os
import pandas as pd
import numpy as np
from itertools import combinations
import pysam
import statistics

from helpers.get_data import get_feature_data

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"
features_all = None

for cohort in os.listdir(directory):
    file = os.path.join(directory, cohort, cohort + "_all.tsv")
    features = pd.read_csv(file, sep="\t", header=0)
    features['IC50 mean'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].mean(axis=1)
    features['IC50 median'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].median(axis=1)
    
    if features_all is None:
        features_all = features
    else:
        features_all = pd.concat([features_all, features])
features_all = features_all.reset_index() 
netmhc_binding = features_all['NetMHC MT IC50 Score'] <= 500 
netmhcpan_binding = features_all['NetMHCpan MT IC50 Score'] <= 500
mhcflurry_binding = features_all['MHCflurry MT IC50 Score'] <= 500

print("%NetMHC <= 500:", netmhc_binding.sum() / features_all.shape[0])
print("%NetMHCpan <= 500:", netmhcpan_binding.sum() / features_all.shape[0])
print("%MHCflurry <= 500:", mhcflurry_binding.sum() / features_all.shape[0])

print("-> filter by mean would kick:", (features_all['IC50 mean'] <= 500).sum() / features_all.shape[0])
print("-> filter by median would kick:", (features_all['IC50 median'] <= 500).sum() / features_all.shape[0])

intersection_all = netmhc_binding & netmhcpan_binding & mhcflurry_binding
union_all = netmhc_binding | netmhcpan_binding | mhcflurry_binding
jaccard_all = intersection_all.sum() / union_all.sum()

print("Jaccard Index of all binding tools:", jaccard_all)

for (binding1, name1), (binding2, name2) in combinations(zip([netmhc_binding, netmhcpan_binding, mhcflurry_binding], ['NetMHC', 'NetMHCpan', 'MHCflurry']), 2):
    intersection = binding1 & binding2
    union = binding1 | binding2
    jaccard = intersection.sum() / union.sum()
    
    print("Jaccard Index of", name1, name2, ":", jaccard)

directory = "/mnt/storage2/users/ahnelll1/master_thesis/data"

num_variants = []
for cohort in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, cohort)):
        continue
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            filename = os.fsdecode(os.path.join(directory, cohort, method, sample))
            variants = [snv for snv in pysam.VariantFile(filename , "r")]
            num_variants.append(len(variants))

print("Median number of variants per Patient:", statistics.median(num_variants))

feature_data = get_feature_data(True, ['AxelMelanomaPhD', 'SomaticAndTreatment'])

# after prefilter
num_epitopes = []
for i in range(np.nanmax(feature_data['sample_id'])):
    num_epitopes.append((feature_data['sample_id'] == i).sum())
    
print("Median number of neoepitopes per Patient:", statistics.median(num_epitopes))