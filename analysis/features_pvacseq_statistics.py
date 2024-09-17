import os
import pandas as pd
from itertools import combinations


directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

for cohort in os.listdir(directory):
    print(cohort)
    file = os.path.join(directory, cohort, cohort + "_all.tsv")
    features = pd.read_csv(file, sep="\t", header=0)
    features['IC50 mean'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].mean(axis=1)
    features['IC50 median'] = features.loc[:, ['NetMHC MT IC50 Score', 'NetMHCpan MT IC50 Score', 'MHCflurry MT IC50 Score']].median(axis=1)
    
    netmhc_binding = features['NetMHC MT IC50 Score'] <= 500 
    netmhcpan_binding = features['NetMHCpan MT IC50 Score'] <= 500
    mhcflurry_binding = features['MHCflurry MT IC50 Score'] <= 500
    
    print("%NetMHC <= 500:", netmhc_binding.sum() / features.shape[0])
    print("%NetMHCpan <= 500:", netmhcpan_binding.sum() / features.shape[0])
    print("%MHCflurry <= 500:", mhcflurry_binding.sum() / features.shape[0])
    
    print("-> filter by mean would kick:", (features['IC50 mean'] <= 500).sum() / features.shape[0])
    print("-> filter by median would kick:", (features['IC50 median'] <= 500).sum() / features.shape[0])
    
    intersection_all = netmhc_binding & netmhcpan_binding & mhcflurry_binding
    union_all = netmhc_binding | netmhcpan_binding | mhcflurry_binding
    jaccard_all = intersection_all.sum() / union_all.sum()
    
    print("Jaccard Index of all binding tools:", jaccard_all)
    
    for (binding1, name1), (binding2, name2) in combinations(zip([netmhc_binding, netmhcpan_binding, mhcflurry_binding], ['NetMHC', 'NetMHCpan', 'MHCflurry']), 2):
        intersection = binding1 & binding2
        union = binding1 | binding2
        jaccard = intersection.sum() / union.sum()
        
        print("Jaccard Index of", name1, name2, ":", jaccard)
    
    