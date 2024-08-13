# Author: Lena Nellinger
# program: cohort_identification_pipeline.py
# purpose: prepares data, runs pVACseq and annotates with neofox for whole cohort_identification_pipeline

from optparse import OptionParser
import sys
import subprocess
import pandas as pd
import os
import json
import time
import shutil
from shutil import which
import logging

from neofox.references.references import ReferenceFolder
from neofox.model.factories import PatientFactory
from neofox.model.factories import NeoepitopeFactory
from neofox.neofox_epitope import NeoFoxEpitope
from neofox.model.conversion import ModelConverter

from data_preparation.annotate_vep_for_pvacseq import annotate_vep_for_pvacseq
from data_preparation.prepare_for_pvacseq import prepare_input_for_pvacseq

def main(argv):
    if len(argv) == 1:
        argv.append("--help")
    usage = "usage: python cohort_identification_pipeline.py -i vcf_file --strelka --regions intersection_all_srr_targets.bed --pairings pairings.tsv --data_dir /mnt/storage1/projects/research/AxelMelanomaPhD"
    desc = "Prepares data and runs pVACseq."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-i", "--input", action="store", dest="input", type="string", help="Input vcf file")
    parser.add_option("--strelka", action="store_true", dest="strelka", default=False, help="Input is produced by strelka")
    parser.add_option("--dragen", action="store_true", dest="dragen", default=False, help="Input is produced by dragen")
    parser.add_option("--data_dir", action="store", dest="data_dir", type="string", help="RNA and HLA data folder")
    parser.add_option("--regions", action="store", dest="regions", type="string",
                      help="Target regions file for filtering", default="/mnt/storage2/users/ahnelll1/master_thesis/data/intersection_srr_inhouse.bed")
    parser.add_option("--pairings", action="store", dest="pairings", type="string",
                      help="Pairings file for RNA pairing", default="/mnt/storage2/users/ahnelll1/master_thesis/data/all_pairings.tsv")
    parser.add_option("--log", action="store_true", dest="log", help="If logging should be on", default=False)
    (options, args) = parser.parse_args()
    

    
    os.environ["MHCFLURRY_DOWNLOADS_DIR"] = "/mnt/storage2/users/ahnelll1/tools/MHCflurry/4/2.0.0"
    os.environ["NEOFOX_MAKEBLASTDB"] = "/mnt/storage2/users/ahnelll1/tools/ncbi-blast-2.10.1+/bin/makeblastdb"
    os.environ["NEOFOX_BLASTP"] = "/mnt/storage2/users/ahnelll1/tools/ncbi-blast-2.10.1+/bin/blastp"
    os.environ["NEOFOX_RSCRIPT"] = which('Rscript')
    os.environ["NEOFOX_NETMHCPAN"] = "/mnt/storage2/users/ahnelll1/tools/netMHCpan-4.1/netMHCpan"
    os.environ["NEOFOX_NETMHC2PAN"] = "/mnt/storage2/users/ahnelll1/tools/netMHCIIpan-4.0/netMHCIIpan"
    os.environ["NEOFOX_MIXMHCPRED"] = "/mnt/storage2/users/ahnelll1/tools/MixMHCpred-2.2/MixMHCpred"
    os.environ["NEOFOX_MIXMHC2PRED"] = "/mnt/storage2/users/ahnelll1/tools/MixMHC2pred-2.0.2.2/MixMHC2pred_unix"
    os.environ["NEOFOX_PRIME"] = "/mnt/storage2/users/ahnelll1/tools/PRIME-2.0/PRIME"
    os.environ["NEOFOX_REFERENCE_FOLDER"] = "/mnt/storage2/users/ahnelll1/tools/neofox-reference"
    
    if not options.strelka and not options.dragen:
        raise ValueError('Input type should be strelka or dragen.')
        
    reference_folder = ReferenceFolder(organism='human')
    hla_database = reference_folder.get_mhc_database()
    
    filename = os.path.basename(options.input)
    f = options.input
    
    output_dir = os.path.join('/mnt/storage2/users/ahnelll1/master_thesis/output_background', os.path.basename(options.data_dir), 'strelka' if options.strelka else 'dragen', filename.replace('.vcf.gz', ''))
    os.makedirs(output_dir, exist_ok=True)
    
    logging.getLogger('asyncio').setLevel(logging.WARNING)
    if options.log:
        logfile = os.path.join(output_dir, "log.txt")
        logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    else:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    
    # checking if it is a file
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

    # prepare input for pvac_seq adn annotate with expression data
    output_ready_filename = output_vep_filename.replace('.vcf.gz', '_ready.vcf.gz')
    if os.path.exists(output_ready_filename):
        logging.info('2. Preparation for pVACseq and RNA Annotation already done')
    else:
        logging.info('2. Preparation for pVACseq and RNA Annotation running...')
        input_type = 'strelka' if options.strelka else 'dragen' if options.dragen else ''
        prepare_input_for_pvacseq(output_vep_filename, input_type, sample_name, normal_name, output_ready_filename, options.pairings, options.data_dir)

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

    # run pVACseq        
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
                              
    # run neofox
    output_dir_neofox = os.path.join(output_dir, 'neofox')
    if os.path.exists(output_dir_neofox):
        logging.info('5. neofox annotation already done')
    else:
        logging.info('5. neofox annotation running...')
        os.makedirs(output_dir_neofox, exist_ok=True)
        pvacseq_filtered = pd.read_csv(os.path.join(output_dir_pvacseq, 'MHC_Class_I', sample_name + '.filtered.tsv'), sep='\t', header=0)
        
        # create a neoepitope candidate using the factory
        # no RNA data given -> only imputed gene expression
        neoepitopes = [NeoepitopeFactory.build_neoepitope(
            mutated_peptide=neoepitope['MT Epitope Seq'],
            wild_type_peptide=neoepitope['WT Epitope Seq'] if neoepitope['Variant Type'] == 'missense' else '',
            patient_identifier="Ptx",
            gene=neoepitope['Gene Name'],
            allele_mhc_i=neoepitope['HLA Allele'],
            organism='human',
            mhc_database=hla_database,
            rna_variant_allele_frequency=None if pd.isna(neoepitope['Tumor RNA VAF']) else neoepitope['Tumor RNA VAF'],
            dna_variant_allele_frequency=None if pd.isna(neoepitope['Tumor DNA VAF']) else neoepitope['Tumor DNA VAF'],
            Chromosome=neoepitope['Chromosome'],
            Start=neoepitope['Start'],
            Stop=neoepitope['Stop'],
            Transcript=neoepitope['Transcript']
        ) for _, neoepitope in pvacseq_filtered.iterrows()]

        patient = PatientFactory.build_patient(
            identifier="Ptx",
            mhc_alleles=hla_alleles,
            mhc2_alleles=[],
            mhc_database=hla_database,
            tumor_type="SKCM"
        )
        annotated_neoepitopes = NeoFoxEpitope(neoepitopes=neoepitopes, patients=[patient], num_cpus=5, log_file_name=logfile if options.log else None).get_annotations()
        annotations_table = ModelConverter.annotated_neoepitopes2epitopes_table(neoepitopes=annotated_neoepitopes, mhc='mhcI')
        annotations_table.to_csv(os.path.join(output_dir_neofox, sample_name + "_neofox_annotations.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main(sys.argv)
