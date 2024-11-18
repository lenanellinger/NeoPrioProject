import pandas as pd
import shutil
import os
import sys
from optparse import OptionParser


def main(argv):
    usage = "usage: python srr_remove_duplicates.py --sample-dir /mnt/storage2/users/ahnelll1/master_thesis/data/srr_vcfs/strelka --metadata /mnt/storage2/users/ahnelll1/master_thesis/data/srr_metadata.tsv"
    desc = "Only keeps last sample per patient if multiple samples for cohort of Axel."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--sample-dir", action="store", dest="sample_dir",
                      help="Directory with sample VCF files")
    parser.add_option("--metadata", action="store", dest="metadata",
                      help="Metadata tsv file.")
    (options, args) = parser.parse_args()

    metadata = pd.read_csv(options.metadata, sep="\t", header=0)
    duplicates = metadata[(metadata["COHORT"] + metadata["PATIENT"]).duplicated('last')]

    print(duplicates)
    for _, d in duplicates.iterrows():
        file = os.path.join(options.sample_dir, d['TUMOR'] + "-" + d['NORMAL'] + '.vcf.gz')
        if os.path.isfile(file):
            print("duplicated: move", file)
            os.makedirs(os.path.join(options.sample_dir, "_no_quality"), exist_ok=True)
            shutil.move(file, os.path.join(options.sample_dir, "_no_quality", d['TUMOR'] + "-" + d['NORMAL'] + '.vcf.gz'))
        else:
            print("duplicated but no file:", file)


if __name__ == "__main__":
    main(sys.argv)