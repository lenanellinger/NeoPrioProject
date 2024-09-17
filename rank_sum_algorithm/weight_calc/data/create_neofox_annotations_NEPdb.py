import pandas as pd
import os
from shutil import which

from neofox.references.references import ReferenceFolder
from neofox.model.factories import PatientFactory
from neofox.model.factories import NeoepitopeFactory
from neofox.neofox_epitope import NeoFoxEpitope
from neofox.model.conversion import ModelConverter

def main():

    os.environ["NEOFOX_MAKEBLASTDB"] = "/mnt/storage2/users/ahnelll1/tools/ncbi-blast-2.10.1+/bin/makeblastdb"
    os.environ["NEOFOX_BLASTP"] = "/mnt/storage2/users/ahnelll1/tools/ncbi-blast-2.10.1+/bin/blastp"
    os.environ["NEOFOX_RSCRIPT"] = which('Rscript')
    os.environ["NEOFOX_NETMHCPAN"] = "/mnt/storage2/users/ahnelll1/tools/netMHCpan-4.1/netMHCpan"
    os.environ["NEOFOX_NETMHC2PAN"] = "/mnt/storage2/users/ahnelll1/tools/netMHCIIpan-4.0/netMHCIIpan"
    os.environ["NEOFOX_MIXMHCPRED"] = "/mnt/storage2/users/ahnelll1/tools/MixMHCpred-2.2/MixMHCpred"
    os.environ["NEOFOX_MIXMHC2PRED"] = "/mnt/storage2/users/ahnelll1/tools/MixMHC2pred-2.0.2.2/MixMHC2pred_unix"
    os.environ["NEOFOX_PRIME"] = "/mnt/storage2/users/ahnelll1/tools/PRIME-2.0/PRIME"
    os.environ["NEOFOX_REFERENCE_FOLDER"] = "/mnt/storage2/users/ahnelll1/tools/neofox-reference"
        
    reference_folder = ReferenceFolder(organism='human')
    hla_database = reference_folder.get_mhc_database()

    samples = pd.read_csv('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/weight_calc/data/NECID_Query.csv')

    # filter length and MHC Class I and multiple identical samples (with different TCR information)
    samples = samples[(samples['antigen_len'] >= 8) & (samples['antigen_len'] <= 14) & (~pd.isna(samples['alleleA'])) & (~samples["alleleA"].str.startswith("DQA", na=True))]
    samples = samples.drop_duplicates(subset=['mut_peptide', 'wt_peptide', 'genesymbol', 'alleleA', 'Tumor Type'], keep='first')

    tumor_type_mapping = {
        'Lung cancer': 'LUAD', 
        'Breast cancer': 'BRCA', 
        'Ovarian cancer': 'OV', 
        'Cervical Carcinoma(HPV+)': 'CESC', 
        'Melanoma': 'SKCM', 
        'SC Lung cancer': 'LUAD',
        'Esophageal adenocarcinoma': 'ESCA', 
        'Uterine carcinoma': 'UCEC', 
        'Endometrial cancer': 'UCEC', 
        'Colorectal cancer': 'READ', 
        'NSC Lung cancer': 'LUSC',
        'Bladder cancer': 'BLCA', 
        'Pancreatic cancer': 'PAAD', 
        'Head and Neck cancer': 'HNSC'
    }

    neoepitopes = [NeoepitopeFactory.build_neoepitope(
        mutated_peptide=neoepitope['mut_peptide'],
        wild_type_peptide=neoepitope['wt_peptide'] if len(neoepitope['wt_peptide']) == len(neoepitope['mut_peptide']) else None,
        patient_identifier=str(neoepitope['id']),
        gene=neoepitope['genesymbol'],
        allele_mhc_i=neoepitope['alleleA'],
        organism='human',
        mhc_database=hla_database,
        response=neoepitope['response']
    ) for _, neoepitope in samples.iterrows()]

    patients = [PatientFactory.build_patient(
        identifier=str(neoepitope['id']),
        mhc_alleles=[neoepitope['alleleA']],
        mhc2_alleles=[],
        mhc_database=hla_database,
        tumor_type=tumor_type_mapping[neoepitope['Tumor Type']] if neoepitope['Tumor Type'] in tumor_type_mapping else None
    ) for _, neoepitope in samples.iterrows()]

    annotated_neoepitopes = NeoFoxEpitope(neoepitopes=neoepitopes, patients=patients, num_cpus=1).get_annotations()
    annotations_table = ModelConverter.annotated_neoepitopes2epitopes_table(neoepitopes=annotated_neoepitopes, mhc='mhcI')
    annotations_table.to_csv(os.path.join("/mnt/storage2/users/ahnelll1/master_thesis/training_data", "NEPdb_neofox_annotations.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()