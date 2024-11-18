import subprocess
import logging
    
def annotate_vep_for_pvacseq(input_file, output_file, vep_command="/mnt/storage2/megSAP/tools/ensembl-vep-release-110.1/vep", threads=5, reference="/mnt/storage2/megSAP/data/genomes/GRCh38.fa", cache="/tmp/local_ngs_data_GRCh38/ensembl-vep-110/", plugins="/mnt/storage2/users/ahgscha1/phd/tools/VEP_plugins/"):
    """
    Annotates VCF file with VEP-data according pVACseq requirements
    
    Keyword arguments:
    input_file -- Input somatic VCF
    output_file -- Output VEP-annotated vcf.gz
    vep_command -- Path to VEP instance (default: /mnt/storage2/megSAP/tools/ensembl-vep-release-110.1/vep)
    threads -- Number of threads to be used (default: 5)
    reference -- Path to reference genome (default: /mnt/storage2/megSAP/data/genomes/GRCh38.fa)
    cache -- Path to VEP cache dir (default: /tmp/local_ngs_data_GRCh38/ensembl-vep-110/)
    plugins -- Path to VEP plugin dir (default: /mnt/storage2/users/ahgscha1/phd/tools/VEP_plugins/)
    """
    vep_command = [
            vep_command,
            "--fork", str(threads),
            "--input_file", input_file,
            "--output_file", output_file.replace(".gz",""),
            "--fasta", reference,
            "--dir_cache", cache,
            "--dir_plugins", plugins,
            "--plugin", "Frameshift",
            "--plugin", "Wildtype",
            "--format", "vcf",
            "--vcf", "--symbol", "--terms", "SO",
            "--tsl", "--biotype", "--hgvs",
            "--offline", "--cache"
        ]

    proc = subprocess.run(vep_command,
        capture_output=True,
        text=True,
        check=False)
    for errline in proc.stderr.splitlines():
        logging.debug("%s: %s", "annotate_vep_for_pvacseq", errline.rstrip())
 
    proc = subprocess.run(["bgzip", output_file.replace(".gz","")],
                                capture_output=True,
                                text=True,
                                check=False)
    for errline in proc.stderr.splitlines():
        logging.debug("%s: %s", "annotate_vep_for_pvacseq", errline.rstrip())

    proc = subprocess.run(["tabix", "-f", "-p", "vcf", output_file],
                                capture_output=True,
                                text=True,
                                check=False)
    for errline in proc.stderr.splitlines():
        logging.debug("%s: %s", "annotate_vep_for_pvacseq", errline.rstrip())
