import os
import pandas as pd
import shutil

"""
Checks if for all commands in commands_quality_filtered.sh a valid output was generated or if an error occurred.
Commands which are not associated with a valid output are written in commands_missing.sh and can be retried.
"""

directory = "/mnt/storage2/users/ahnelll1/master_thesis/output"

with open('/mnt/storage2/users/ahnelll1/master_thesis/pipeline_commands/commands_quality_filtered.sh', 'r') as fi:
    with open('/mnt/storage2/users/ahnelll1/master_thesis/pipeline_commands/commands_missing.sh', 'w') as fo:
        lines = fi.readlines()
        for cohort in os.listdir(directory):
            for method in os.listdir(os.path.join(directory, cohort)):
                if os.path.isfile(os.path.join(directory, cohort, method)):
                    continue
                out = []
                for sample in os.listdir(os.path.join(directory, cohort, method)):
                    if sample.startswith("_no"):
                        continue
                    file = os.path.join(directory, cohort, method, sample)
                        
                    command = [l for l in lines if sample in l]
                    if len(command) != 1:
                        print("Not exactly one command found:", sample)
                        continue
                    else:
                        command = command[0]
                    if os.path.isfile(os.path.join(file, "pVACseq", "MHC_Class_I", sample.split("-")[0] + ".filtered.tsv")):
                        tsv = pd.read_csv(os.path.join(file, "pVACseq", "MHC_Class_I", sample.split("-")[0] + ".filtered.tsv"), sep="\t", header=0)
                        if not os.path.isdir(os.path.join(file, "neofox")):
                            # no neofox folder
                            fo.write(command)
                            out.append(sample)
                        elif len([name for name in os.listdir(os.path.join(file, "neofox"))]) == 0 and tsv.shape[0] > 0:
                            # neoepitopes in pVACseq output, but no neofox output file
                            fo.write(command)
                            out.append(sample)
                            shutil.rmtree(os.path.join(directory, cohort, method, sample, 'neofox'))
                    else:
                        # No pVACseq output
                        if os.path.isdir(os.path.join(file, "neofox")):
                            shutil.rmtree(os.path.join(directory, cohort, method, sample, 'neofox'))
                        else:
                            print("No pVAC and no neofox:", sample)
                        fo.write(command)
                        out.append(sample)
                
                
                print(cohort, method, out)
                print(len(out))