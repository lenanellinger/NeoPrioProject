# Author: Lena Nellinger

from optparse import OptionParser
import sys
import pandas as pd
import os
from shutil import which
import logging

from neofox.references.references import ReferenceFolder
from neofox.model.factories import PatientFactory
from neofox.model.factories import NeoepitopeFactory
from neofox.neofox_epitope import NeoFoxEpitope
from neofox.model.conversion import ModelConverter


def main(argv):
    if len(argv) == 1:
        argv.append("--help")
    usage = "usage: python run_neofox.py -s <sample_nmae> -n <normal_name> --strelka --data_dir /mnt/storage1/projects/research/AxelMelanomaPhD"
    desc = "Annotates with neofox."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-s", "--sample", action="store", dest="sample", type="string", help="Sample Name")
    parser.add_option("-n", "--normal", action="store", dest="normal", type="string", help="Normal Name")
    parser.add_option("--strelka", action="store_true", dest="strelka", default=False,
                      help="Input is produced by strelka")
    parser.add_option("--dragen", action="store_true", dest="dragen", default=False, help="Input is produced by dragen")
    parser.add_option("--data_dir", action="store", dest="data_dir", type="string", help="RNA and HLA data folder")
    parser.add_option("--log", action="store_true", dest="log", help="If logging should be on", default=False)
    (options, args) = parser.parse_args()

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

    output_dir = os.path.join('/mnt/storage2/users/ahnelll1/master_thesis/output', os.path.basename(options.data_dir),
                              'strelka' if options.strelka else 'dragen', options.sample + '-' + options.normal)
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
    logger.info('Processing file ' + options.sample + '-' + options.normal)

    # prepare data
    sample_name = options.sample
    normal_name = options.normal

    # run neofox
    output_dir_neofox = os.path.join(output_dir, 'neofox')
    hla_file = os.path.join(options.data_dir, "Sample_" + normal_name, normal_name + '_hla_genotyper.tsv')
    if not os.path.isfile(hla_file):
        logging.warning("HLA file does not exists")
        exit(0)
    alleles_df = pd.read_csv(hla_file, sep='\t', header=0)
    hla_alleles = [row['a1'] for _, row in alleles_df.iterrows()] + [row['a2'] for _, row in alleles_df.iterrows()]
    if os.path.exists(output_dir_neofox):
        logging.info('5. neofox annotation already done')
    else:
        logging.info('5. neofox annotation running...')
        os.makedirs(output_dir_neofox, exist_ok=True)
        pvacseq_filtered = pd.read_csv(os.path.join(output_dir, 'pVACseq', 'MHC_Class_I', sample_name + '.filtered.tsv'),
                                       sep='\t', header=0)

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
        annotated_neoepitopes = NeoFoxEpitope(neoepitopes=neoepitopes, patients=[patient], num_cpus=5,
                                              log_file_name=logfile if options.log else None).get_annotations()
        annotations_table = ModelConverter.annotated_neoepitopes2epitopes_table(neoepitopes=annotated_neoepitopes,
                                                                                mhc='mhcI')
        annotations_table.to_csv(os.path.join(output_dir_neofox, sample_name + "_neofox_annotations.tsv"), sep="\t",
                                 index=False)


if __name__ == "__main__":
    main(sys.argv)
