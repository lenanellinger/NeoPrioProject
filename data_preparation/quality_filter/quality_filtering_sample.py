import os
import shutil
import xml.etree.ElementTree as ET

# define directory names specific for your interest
directory = "/mnt/storage2/users/ahnelll1/master_thesis/output/AxelMelanomaPhD/strelka"
directory_samples = "/mnt/storage1/projects/research/AxelMelanomaPhD"

def is_quality_parameter_sufficient(qcml_file, parameter, threhsold):
    tree = ET.parse(qcml_file)
    root = tree.getroot()
    sufficient = True

    for item in root.findall(
            './{http://www.prime-xs.eu/ms/qcml}runQuality/{http://www.prime-xs.eu/ms/qcml}qualityParameter'):
        if item.get('name') == parameter and float(item.get('value')) < threhsold:
            # if you want to move output files which do not meet the quality criteria, uncomment the following line
            # shutil.move(file, os.path.join(directory, "_no_quality", filename))
            print(filename)
            sufficient = False
            break

    return sufficient

for filename in os.listdir(directory):
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

    if not is_quality_parameter_sufficient(qcml_file, "target region read depth", 60.):
        continue
        
    # Coverage Normal
    qcml_file = os.path.join(directory_samples, "Sample_" + normal_id, normal_id + "_stats_map.qcML")
    if not os.path.isfile(qcml_file):
        print("No folder for normal " + normal_id)
        continue

    if not is_quality_parameter_sufficient(qcml_file, "target region read depth", 20.):
        continue
    
    # Correlation Tumor - Normal
    qcml_file = os.path.join(directory_samples, "Somatic_" + filename, filename + "_stats_som.qcML")
    if not os.path.isfile(qcml_file):
        print("No folder for somatic " + filename)
        continue

    if is_quality_parameter_sufficient(qcml_file, "	sample correlation", 0.8):
        continue
        
    # Multiple Samples
    if not tumor_id.endswith("_01") or not normal_id.endswith("_01"):
        print(tumor_id, normal_id)