import os
import shutil
import xml.etree.ElementTree as ET

# for both cohorts and each for strelka und dragen
directory = "/mnt/storage2/users/ahnelll1/master_thesis/output_background/AxelMelanomaPhD/strelka"
directory_samples = "/mnt/storage1/projects/research/AxelMelanomaPhD"

for filename in os.listdir(directory):
    stop = False
    
    if filename.startswith("_no"):
        continue

    file = os.path.join(directory, filename)
    
    tumor_id = filename.split("-")[0]
    normal_id = filename.split("-")[1]
    
    # Coverage Tumor
    qcml_file = os.path.join(directory_samples, "Sample_" + tumor_id, tumor_id + "_stats_map.qcML")
    if not os.path.isfile(qcml_file):
        print("No folder for tumor " + tumor_id)
        continue
    tree = ET.parse(qcml_file)
    root = tree.getroot()
            
    for item in root.findall('./{http://www.prime-xs.eu/ms/qcml}runQuality/{http://www.prime-xs.eu/ms/qcml}qualityParameter'):
        if item.get('name') == "target region read depth" and float(item.get('value')) < 60.:
            #shutil.move(file, os.path.join(directory, "_no_quality", filename))
            print(filename)
            stop = True
            break
    if stop:
        continue
        
    # Coverage Normal
    qcml_file = os.path.join(directory_samples, "Sample_" + normal_id, normal_id + "_stats_map.qcML")
    if not os.path.isfile(qcml_file):
        print("No folder for normal " + normal_id)
        continue
    tree = ET.parse(qcml_file)
    root = tree.getroot()
            
    for item in root.findall('./{http://www.prime-xs.eu/ms/qcml}runQuality/{http://www.prime-xs.eu/ms/qcml}qualityParameter'):
        if item.get('name') == "target region read depth" and float(item.get('value')) < 20.:
            #shutil.move(file, os.path.join(directory, "_no_quality", filename))
            print(filename)
            stop = True
            break
    if stop:
        continue
    
    # Correlation
    qcml_file = os.path.join(directory_samples, "Somatic_" + filename, filename + "_stats_som.qcML")
    if not os.path.isfile(qcml_file):
        print("No folder for somatic " + filename)
        continue
    tree = ET.parse(qcml_file)
    root = tree.getroot()
            
    for item in root.findall('./{http://www.prime-xs.eu/ms/qcml}runQuality/{http://www.prime-xs.eu/ms/qcml}qualityParameter'):
        if item.get('name') == "	sample correlation" and float(item.get('value')) < 0.8:
            #shutil.move(file, os.path.join(directory, "_no_quality", filename))
            print(filename)
            stop = True
            break
    if stop:
        continue
        
    # Multiple Samples
    if not tumor_id.endswith("_01") or not normal_id.endswith("_01"):
        print(tumor_id, normal_id)