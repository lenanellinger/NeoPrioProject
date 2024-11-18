import sys
import shutil
import os
import xml.etree.ElementTree as ET
from optparse import OptionParser

def is_quality_parameter_sufficient(qcml_file, parameter, threhsold, dir, filename):
    tree = ET.parse(qcml_file)
    root = tree.getroot()
    sufficient = True

    for item in root.findall(
            './{http://www.prime-xs.eu/ms/qcml}runQuality/{http://www.prime-xs.eu/ms/qcml}qualityParameter'):
        if item.get('name') == parameter and float(item.get('value')) < threhsold:
            os.makedirs(os.path.join(dir, "_no_quality"), exist_ok=True)
            shutil.move(os.path.join(dir, filename), os.path.join(dir, "_no_quality", filename))
            print(filename)
            sufficient = False
            break

    return sufficient


def main(argv):
    usage = "usage: python quality_filtering_sample.py --sample-dir /mnt/storage2/users/ahnelll1/master_thesis/data/srr_vcfs/strelka --seq-dir /mnt/storage1/projects/research/AxelMelanomaPhD"
    desc = "Quality filters on sample level."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--sample-dir", action="store", dest="sample_dir",
                      help="Directory with sample VCF files")
    parser.add_option("--seq-dir", action="store", dest="seq_dir", help="Directory with sequencing data for each sample")
    (options, args) = parser.parse_args()

    for filename in os.listdir(options.sample_dir):
        if not os.path.isfile(os.path.join(options.sample_dir, filename)):
            continue

        tumor_id = filename.split("-")[0]
        normal_id = filename.split("-")[1].replace(".vcf.gz", "")

        # Coverage Tumor
        qcml_file = os.path.join(options.seq_dir, "Sample_" + tumor_id, tumor_id + "_stats_map.qcML")
        if not os.path.isfile(qcml_file):
            print("No folder for tumor " + tumor_id)
            continue

        if not is_quality_parameter_sufficient(qcml_file, "target region read depth", 60., options.sample_dir, filename):
            continue

        # Coverage Normal
        qcml_file = os.path.join(options.seq_dir, "Sample_" + normal_id, normal_id + "_stats_map.qcML")
        if not os.path.isfile(qcml_file):
            print("No folder for normal " + normal_id)
            continue

        if not is_quality_parameter_sufficient(qcml_file, "target region read depth", 20., options.sample_dir, filename):
            continue

        # Correlation Tumor - Normal
        qcml_file = os.path.join(options.seq_dir, "Somatic_" + filename, filename + "_stats_som.qcML")
        if not os.path.isfile(qcml_file):
            print("No folder for somatic " + filename)
            continue

        if is_quality_parameter_sufficient(qcml_file, "	sample correlation", 0.8, options.sample_dir, filename):
            continue

        # Multiple Samples
        if not tumor_id.endswith("_01") or not normal_id.endswith("_01"):
            os.makedirs(os.path.join(options.sample_dir, "_no_quality"), exist_ok=True)
            shutil.move(os.path.join(options.sample_dir, filename), os.path.join(options.sample_dir, "_no_quality", filename))
            print(filename)




if __name__ == "__main__":
    main(sys.argv)