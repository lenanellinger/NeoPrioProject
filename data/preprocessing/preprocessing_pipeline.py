# Author: Lena Nellinger

from optparse import OptionParser
import sys
import subprocess
import os
import logging

from annotate_vep_for_pvacseq import annotate_vep_for_pvacseq
from prepare_for_pvacseq import prepare_input_for_pvacseq


def main(argv):
    if len(argv) == 1:
        argv.append("--help")
    usage = "usage: python preprocessing_pipeline.py -i vcf_file --strelka --regions intersection_all_srr_targets.bed --pairings pairings.tsv --data_dir /mnt/storage1/projects/research/AxelMelanomaPhD"
    desc = "Preprocesses data for pVACseq. Quality filtering not part of the pipeline."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-i", "--input", action="store", dest="input", type="string", help="Input vcf file")
    parser.add_option("--strelka", action="store_true", dest="strelka", default=False,
                      help="Input is produced by strelka")
    parser.add_option("--dragen", action="store_true", dest="dragen", default=False, help="Input is produced by dragen")
    parser.add_option("--data_dir", action="store", dest="data_dir", type="string", help="RNA and HLA data folder")
    parser.add_option("--regions", action="store", dest="regions", type="string",
                      help="Target regions file for filtering",
                      default="/mnt/storage2/users/ahnelll1/master_thesis/data/intersection_srr_inhouse.bed")
    parser.add_option("--pairings", action="store", dest="pairings", type="string",
                      help="Pairings file for RNA pairing",
                      default="/mnt/storage2/users/ahnelll1/master_thesis/data/all_pairings.tsv")
    parser.add_option("--log", action="store_true", dest="log", help="If logging should be on", default=False)
    (options, args) = parser.parse_args()

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

    # annotate with VEP
    output_vep_filename = os.path.join(output_dir, filename.replace('.vcf.gz', '_vep.vcf.gz'))
    if os.path.exists(output_vep_filename):
        logging.info('1. VEP already annotated')
    else:
        logging.info('1. VEP Annotation running...')
        annotate_vep_for_pvacseq(f, output_vep_filename)

    # prepare input for pVACseq and annotate with expression data
    output_ready_filename = output_vep_filename.replace('.vcf.gz', '_ready.vcf.gz')
    if os.path.exists(output_ready_filename):
        logging.info('2. Preparation for pVACseq and RNA Annotation already done')
    else:
        logging.info('2. Preparation for pVACseq and RNA Annotation running...')
        input_type = 'strelka' if options.strelka else 'dragen' if options.dragen else ''
        prepare_input_for_pvacseq(output_vep_filename, input_type, sample_name, normal_name, output_ready_filename,
                                  options.pairings, options.data_dir)

    # filter target regions
    output_filter_filename = output_ready_filename.replace('.vcf.gz', '_filtered.vcf.gz')
    if os.path.exists(output_filter_filename):
        logging.info('3. Filter target regions already done')
    else:
        logging.info('3. Filter target regions running...')
        output_filter_file = open(output_filter_filename, 'wb')  # output file
        zcat_proc = subprocess.Popen(['zcat', output_ready_filename], stdout=subprocess.PIPE, text=True)
        filter_reg_proc = subprocess.Popen(['VcfFilter', '-reg', options.regions], stdin=zcat_proc.stdout,
                                           stdout=subprocess.PIPE, text=True, shell=True)
        filter_dt_proc = subprocess.Popen(['VcfFilter -filter_exclude depth-tum'], stdin=filter_reg_proc.stdout,
                                          stdout=subprocess.PIPE, text=True, shell=True)
        filter_dn_proc = subprocess.Popen(['VcfFilter -filter_exclude depth-nor'], stdin=filter_dt_proc.stdout,
                                          stdout=subprocess.PIPE, text=True, shell=True)
        filter_ft_proc = subprocess.Popen(['VcfFilter -filter_exclude freq-tum'], stdin=filter_dn_proc.stdout,
                                          stdout=subprocess.PIPE, text=True, shell=True)
        filter_fn_proc = subprocess.Popen(['VcfFilter -filter_exclude freq-nor | gzip -c'], shell=True,
                                          stdin=filter_ft_proc.stdout, stdout=output_filter_file)


if __name__ == "__main__":
    main(sys.argv)
