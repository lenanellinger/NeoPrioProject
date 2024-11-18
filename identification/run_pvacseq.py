# Author: Lena Nellinger

from optparse import OptionParser
import sys
import subprocess
import pandas as pd
import os
import logging


def main(argv):
    if len(argv) == 1:
        argv.append("--help")
    usage = "usage: python run_pvacseq.py -i vcf_file --strelka --regions intersection_all_srr_targets.bed --pairings pairings.tsv --data_dir /mnt/storage1/projects/research/AxelMelanomaPhD"
    desc = "Runs pVACseq on preprocessed data."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-i", "--input", action="store", dest="input", type="string", help="Input vcf file")
    parser.add_option("--strelka", action="store_true", dest="strelka", default=False,
                      help="Input is produced by strelka")
    parser.add_option("--dragen", action="store_true", dest="dragen", default=False, help="Input is produced by dragen")
    parser.add_option("--data_dir", action="store", dest="data_dir", type="string", help="RNA and HLA data folder")
    parser.add_option("--log", action="store_true", dest="log", help="If logging should be on", default=False)
    (options, args) = parser.parse_args()

    os.environ["MHCFLURRY_DOWNLOADS_DIR"] = "/mnt/storage2/users/ahnelll1/tools/MHCflurry/4/2.0.0"

    if not options.strelka and not options.dragen:
        raise ValueError('Input type should be strelka or dragen.')

    filename = os.path.basename(options.input)
    f = options.input

    output_dir = os.path.join('/mnt/storage2/users/ahnelll1/master_thesis/output', os.path.basename(options.data_dir),
                              'strelka' if options.strelka else 'dragen', filename.replace('.vcf.gz', ''))
    os.makedirs(output_dir, exist_ok=True)

    logging.getLogger('asyncio').setLevel(logging.WARNING)
    if options.log:
        logfile = os.path.join(output_dir, "log.txt")
        logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",
                            format="%(asctime)-15s %(levelname)-8s %(message)s")
    else:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)

    # checking if input file is a file
    if not os.path.isfile(f) or not filename.endswith('.vcf.gz'):
        logger.info('Input object ' + filename + ' could not be parsed. Should be .vcf.gz')
        exit(1)
    logger.info('Processing file ' + filename)

    # prepare data
    sample_name = filename.split('-')[0]
    normal_name = filename.split('-')[1].replace('.vcf.gz', '')

    # run pVACseq
    output_filter_filename = os.path.join(output_dir, filename.replace('.vcf.gz', '_vep_ready_filtered.vcf.gz'))
    output_dir_pvacseq = os.path.join(output_dir, 'pVACseq')
    hla_file = os.path.join(options.data_dir, "Sample_" + normal_name, normal_name + '_hla_genotyper.tsv')
    if not os.path.isfile(hla_file):
        logging.warning("HLA file does not exists")
        exit(0)
    alleles_df = pd.read_csv(hla_file, sep='\t', header=0)
    hla_alleles = [row['a1'] for _, row in alleles_df.iterrows()] + [row['a2'] for _, row in alleles_df.iterrows()]
    if os.path.exists(output_dir_pvacseq):
        logging.info('4. pVACseq already done')
    else:
        logging.info('4. pVACseq running...')
        os.makedirs(output_dir_pvacseq, exist_ok=True)

        proc = subprocess.run(['pvacseq', 'run',
                               output_filter_filename, sample_name, ','.join(hla_alleles),
                               'MHCflurry', 'NetMHC', 'NetMHCpan',
                               output_dir_pvacseq,
                               '-e1', '9',
                               '--normal-sample-name', normal_name,
                               '--netmhc-stab',
                               '-m', 'lowest',
                               '-a', 'sample_name',
                               '--iedb-install-directory', '/mnt/storage2/users/ahgscha1/phd/tools/IEDB_MHC',
                               '-t', '5'],
                              capture_output=True,
                              text=True,
                              check=False)
        for errline in proc.stderr.splitlines():
            logging.debug("%s: %s", "pvacseq", errline.rstrip())



if __name__ == "__main__":
    main(sys.argv)
