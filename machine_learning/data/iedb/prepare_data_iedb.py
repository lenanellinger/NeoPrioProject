import pandas as pd
import numpy as np
import os

assays = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), "tcell_table_export_1725363429.tsv"), sep="\t", header=0)

print(assays.head)

filtered_assays = assays[(assays["Related Object - Epitope Relation"] != "unspecified") & (assays["Related Object - Species"] == "Homo sapiens")]

print(filtered_assays.groupby(["Epitope - CEDAR IRI"]).size().max())
