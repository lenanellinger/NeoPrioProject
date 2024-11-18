import os
import pandas as pd
import statistics

all_alleles = []
netmhc = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:02", "HLA-A*02:03", "HLA-A*02:05", "HLA-A*02:06", "HLA-A*02:07", "HLA-A*02:11", "HLA-A*02:12", "HLA-A*02:16", "HLA-A*02:17", "HLA-A*02:19", "HLA-A*02:50", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*24:03", "HLA-A*25:01", "HLA-A*26:01", "HLA-A*26:02", "HLA-A*26:03", "HLA-A*29:02", "HLA-A*30:01", "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*32:07", "HLA-A*32:15", "HLA-A*33:01", "HLA-A*66:01", "HLA-A*68:01", "HLA-A*68:02", "HLA-A*68:23", "HLA-A*69:01", "HLA-A*80:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*08:02", "HLA-B*08:03", "HLA-B*14:02", "HLA-B*15:01", "HLA-B*15:02", "HLA-B*15:03", "HLA-B*15:09", "HLA-B*15:17", "HLA-B*18:01", "HLA-B*27:05", "HLA-B*27:20", "HLA-B*35:01", "HLA-B*35:03", "HLA-B*38:01", "HLA-B*39:01", "HLA-B*40:01", "HLA-B*40:02", "HLA-B*40:13", "HLA-B*42:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*45:01", "HLA-B*46:01", "HLA-B*48:01", "HLA-B*51:01", "HLA-B*53:01", "HLA-B*54:01", "HLA-B*57:01", "HLA-B*58:01", "HLA-B*58:02", "HLA-B*73:01", "HLA-B*83:01", "HLA-C*03:03", "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*07:01", "HLA-C*07:02", "HLA-C*08:02", "HLA-C*12:03", "HLA-C*14:02", "HLA-C*15:02", "HLA-E*01:01"]
print("Number of HLA alleles in NetMHC:", len(netmhc))
print()

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"
data_dir = {"AxelMelanomaPhD": "/mnt/storage1/projects/research/AxelMelanomaPhD", "SomaticAndTreatment": "/mnt/storage2/projects/diagnostic/SomaticAndTreatment"}

num = []

for cohort in os.listdir(directory):
    for method in os.listdir(os.path.join(directory, cohort)):
        if os.path.isfile(os.path.join(directory, cohort, method)):
            continue
        for sample in os.listdir(os.path.join(directory, cohort, method)):
            if sample.startswith("_no"):
                continue
            normal_name = sample.split('-')[1]
            hla_file = os.path.join(data_dir[cohort], "Sample_" + normal_name, normal_name + '_hla_genotyper.tsv')
            if not os.path.isfile(hla_file):
                hla_file = os.path.join("/mnt/storage2/projects/diagnostic/Exome_Diagnostik", "Sample_" + normal_name, normal_name + '_hla_genotyper.tsv')
            alleles_df = pd.read_csv(hla_file, sep='\t', header=0)
            hla_alleles = [row['a1'] for _, row in alleles_df.iterrows()] + [row['a2'] for _, row in alleles_df.iterrows()]
            all_alleles.extend(hla_alleles)
            num.append(len(hla_alleles))
    
dict_al = dict((x, all_alleles.count(x)) for x in sorted(set(all_alleles)))
no_netmhc_count = 0
all_count = 0
for key in dict_al:
    all_count += dict_al[key]
    if key not in netmhc:
        print("NO NETMCHC", key, ":", dict_al[key])
        no_netmhc_count += dict_al[key]
    else:
        print(key, ":", dict_al[key])

print()
print(f"Not covered with NetMHC: {no_netmhc_count}/{all_count}")
print(f"HLA alleles per patient median: {statistics.median(num)}")
    
