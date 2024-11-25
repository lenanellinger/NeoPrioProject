# NeoPrioProject

## Overview
NeoPrioProject is a Python-based project designed to prioritize neoepitopes using a rank sum algorithm based on neofox features.

## Features
- Identifies neoepitope candidates based on binding affinity utilizing pVACseq
- Prioritizes neoepitopes based on various neofox features
- Extracts feature importance based on machine learning classifier
- Utilizes rank sum algorithm and relative quantile-based quality score for scoring
- Analyzes features as well as identification and prioritization outputs

## Installation
To install the necessary dependencies, run:
```bash
pip install -r requirements.txt
```

## Usage
To identify neoepitope candidates for a patient with pVACseq, run the following command: 
```bash
python identification/run_pvacseq.py -i <path_to_vcf_file> --strelka --data_dir <path_to_directory_containing_all_sequencing_data>
```

To annotate neoepitope candidates with neofox, run:
```bash
python prioritization/run_neofox.py -i <vcf_file> --strelka --data_dir <path_to_directory_containing_all_sequencing_data>
```

To prioritize neoepitopes within a patient, run the following command:
```bash
python prioritization/neoantigen_prioritization_rank_sum.py --pvacseq-filtered-output <path_to_pvacseq_filtered_output> --neofox-output <path_to_neofox_output>
```


## Prerequisites
- Ensure that VcfFilter is in your system's PATH