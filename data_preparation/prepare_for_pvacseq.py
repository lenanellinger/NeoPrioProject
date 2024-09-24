import pysam
import pysamstats
import pandas as pd
import os
import logging
variant_types = ["splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
        "missense_variant", "protein_altering_variant", "splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
        "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "NMD_transcript_variant"]

def prepare_input_for_pvacseq(input_file, input_type, sample_name, normal_name, output_file, pairings_tsv, rna_dir):
    """
    Prepares Input for pVACseq and annotates RNA data
    
    Keyword arguments:
    input_file -- input somatic VCF file
    input_type -- strelka or dragen
    sample_name -- name of the sample
    normal_name -- name of normal
    output_file -- VCF file ready for pVACseq
    pairings_tsv -- TSV that pairs sample name, normal name and rna name
    rna_dir -- Path to directory with RNA data
    """
    if input_type not in ['strelka', 'dragen']:
        raise ValueError('Input type should be strelka or dragen')
       
    rna_bam, expr_file = get_RNA(pairings_tsv, sample_name, normal_name, rna_dir)
    
    old_vcf = pysam.VariantFile(input_file , "r")
    add_header(old_vcf, input_type, rna_bam, expr_file)
    new_vcf = pysam.VariantFile(output_file, "w", header=old_vcf.header)
    
    for snv in old_vcf:
        for sample in snv.samples:
            if input_type == 'strelka':
                convert_strelka_input(snv, sample)
            annotate_RNA_data(snv, sample, rna_bam, expr_file)
        new_vcf.write(snv)

    old_vcf.close()
    new_vcf.close()

def get_RNA(pairings_tsv, sample_name, normal_name, rna_dir):
    """
    Return RNA bam and expression file
    
    Keyword arguments:
    pairings_tsv -- TSV that pairs sample name, normal name and rna name
    sample_name -- name of the sample
    normal_name -- name of normal
    rna_dir -- Path to directory with RNA data
    
    Returns: 
    rna_bam -- filename for RNA bam file
    expr_file -- filename for RNA expression file
    """
    rna_bam = ""
    expr_file = ""
    
    # check for RNA
    rna_pairings = pd.read_csv(pairings_tsv, sep='\t', header=0)
    rna_pairing_sample = rna_pairings.loc[
        (rna_pairings['TUMOR'] == sample_name) & (rna_pairings['NORMAL'] == normal_name)]
    if rna_pairing_sample.shape[0] != 1:
        logging.warning(f'For sample {sample_name} not correct amount of rna pairings found: Expression data will not be annotated.')
    else:
        rna_name = rna_pairing_sample.iloc[0]['RNA']
        if not pd.isnull(rna_name):
            logging.info(f'RNA data found: {rna_name}')
            rna_folder = os.path.join(rna_dir, 'Sample_' + rna_name)
            rna_bam = os.path.join(rna_folder, rna_name + '.bam')
            expr_file = os.path.join(rna_folder, rna_name + '_counts.tsv')
    return rna_bam, expr_file
    
def add_header(vcf, input_type, rna_bam, expr_file):
    """
    Adds header format fields
    
    Keyword arguments:
    vcf -- VCF to be added
    input_type -- strelka or dragen
    rna_bam -- filename for RNA bam file
    expr_file -- filename for RNA expression file    
    """
    if input_type == 'strelka':
        vcf.header.formats.add("GT", 1, "Float", "Genotype")
        vcf.header.formats.add("AF", 1, "Float", "Variant Allele Frequency.")
        vcf.header.formats.add("AD", 1, "Integer", "Alternate depth of the SNV.")
            
    if rna_bam != "":
        vcf.header.formats.add("RDP", 1, "Integer", "RNA total read depth")
        vcf.header.formats.add("RAD", 1, "Integer", "RNA alternate allele read depth")
        vcf.header.formats.add("RAF", 1, "Float", "RNA variant allele frequency")
  
    if expr_file != "":
        vcf.header.formats.add("GX", ".", "String", "Gene Expression")


def convert_strelka_input(snv, sample):
    """
    Add Format fields for AF and AD for one sample
    
    Keyword arguments:
    snv -- SNV including many samples
    sample -- the sample to change
    """
    #Add genotype dummy
    snv.samples[sample]["GT"] = (0,1)
    
    #Add tumor/normal AF and AD
    snv_fields = {"AU", "CU", "GU", "TU"}
    indel_fields = {"TIR", "TAR"}
    if snv_fields.issubset( set(snv.samples[sample].keys()) ) :
        (A1,A2) = snv.samples[sample]["AU"]
        (C1,C2) = snv.samples[sample]["CU"]
        (G1,G2) = snv.samples[sample]["GU"]
        (T1, T2) = snv.samples[sample]["TU"]
        counts = {"A" : A1, "C": C1, "G": G1, "T": T1}
        depth =  sum(counts.values())

        obs = snv.alts[0]

        #variant allele frequency
        af = counts[obs] / depth
        #alternate depth
        ad = counts[obs]
        snv.samples[sample]["AF"] = af
        snv.samples[sample]["AD"] = ad
    elif indel_fields.issubset( set(snv.samples[sample].keys()) ):
        #TIR: Reads strongly supporting indel allele for tiers 1,2
        TIR = snv.samples[sample]["TIR"][0]
        #TAR: Reads strongly supporting alternate allele for tiers 1,2
        TAR = snv.samples[sample]["TAR"][0]

        af = 0.
        if TIR+TAR != 0:
            af = TIR/(TIR+TAR)
        snv.samples[sample]["AF"] = af
        snv.samples[sample]["AD"] = TIR
    else:
        logging.error("Could not parse STRELKA variant:", snv)
        exit(1)

def annotate_RNA_data(snv, sample, rna_bam, expr_file):    
    """
    Adds format fields for RNA data, if existing
    
    Keyword arguments:
    snv -- SNV including many samples
    sample -- the sample to change
    rna_bam -- filename for RNA bam file
    expr_file -- filename for RNA expression file
    """
    #Annoate read depth from RNA BAM
    if rna_bam != "":
        rnafile = pysam.AlignmentFile(rna_bam, "rb")
        if not is_exonic(snv.info["CSQ"]):
            snv.samples[sample]["RDP"] = 0
            snv.samples[sample]["RAD"] = 0
            snv.samples[sample]["RAF"] = 0
        else:
            a = pysamstats.stat_variation(rnafile, fafile="/mnt/storage2/megSAP/data/genomes/GRCh38.fa", chrom=snv.chrom, start=snv.start, end=snv.stop, truncate=True, pad=True)
            for rec in a:
                if(len(snv.ref[0]) == 1 and len(snv.alts[0]) == 1):
                    rna_dp = rec["matches"] + rec["mismatches"]
                    rna_alt_dp = rec[snv.alts[0]]
                    rna_af = 0.
                    if rna_dp > 0: 
                        rna_af = rna_alt_dp/rna_dp
                    snv.samples[sample]["RDP"] = rna_dp
                    snv.samples[sample]["RAD"] = rna_alt_dp
                    snv.samples[sample]["RAF"] = rna_af
                else:
                    rna_dp = rec["matches"] + rec["mismatches"]
                    indel_count = rec["insertions"] +  rec["deletions"]
                    rna_af = 0
                    if rna_dp > 0:
                        rna_af = indel_count / rna_dp

                    snv.samples[sample]["RDP"] = rna_dp
                    snv.samples[sample]["RAD"] = indel_count
                    snv.samples[sample]["RAF"] = rna_af

    #Annotate gene expression data
    if expr_file != "":
        expressionfile = pd.read_csv(expr_file, delimiter="\t", index_col=0)
        vep_info = snv.info["CSQ"]

        #annotate expression of first gene in list (pvactools can only parse one value per SNV)
        gene_expression = {}
        for vep in vep_info:
            parts = vep.split("|")
            gene_id = parts[4]
            if not gene_id in expressionfile.index:
                continue
            if gene_id in gene_expression.keys():
                continue
            gene_expression[gene_id] = gene_id + "|" +str(expressionfile.loc[gene_id]["tpm"].item())


        if len(gene_expression) > 0:
            snv.samples[sample]["GX"] = ",".join(gene_expression.values())
            #snv.samples[sample]["GX"] = list( gene_expression.values() )[0] #Take only one, arbitrary first value. See https://github.com/griffithlab/pVACtools/issues/991
    
def is_exonic(vep_info):
    global variant_types
    is_exonic = False

    for vep in vep_info:
        parts = vep.split("|")
        if parts[1] in variant_types:
            is_exonic = True
    return is_exonic