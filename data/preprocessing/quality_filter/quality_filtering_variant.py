import os
import pandas as pd
import sys
from optparse import OptionParser

def main(argv):
    usage = "usage: python quality_filtering_variant.py --pvacseq-dir /mnt/storage2/users/ahnelll1/master_thesis/output/SomaticAndTreatment/dragen"
    desc = "Checks if the variant occurs on at least 3 reads. Other variant specific quality measures are performed by pVACseq."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--pvacseq-dir", action="store", dest="pvacseq_dir",
                      help="Directory with output files of patients by pVACseq")
    (options, args) = parser.parse_args()

    for filename in os.listdir(options.pvacseq_dir):
        if filename == "_no_quality":
            continue
        file = os.path.join(options.pvacseq_dir, filename)

        # minimum 3 read count (= tumor DNA VAF * Depth)
        pvac_file = os.path.join(file, "pVACseq", "MHC_Class_I", filename.split("-")[0] + ".filtered.tsv")
        if not os.path.isfile(pvac_file):
            print(f"No pVACseq file for {filename}")
            continue
        pvacseq = pd.read_csv(pvac_file, sep="\t", header=0)
        for _, variant in pvacseq.iterrows():
            read_count = variant['Tumor DNA VAF'] * variant['Tumor DNA Depth']
            if read_count < 3:
                print(f"Read Count not sufficient for sample {filename} in variant {variant['MT Epitope Seq']}")


if __name__ == "__main__":
    main(sys.argv)