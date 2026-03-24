# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:51:47 2023

@author: copeland
"""


import genomicsTools, os, logging, argparse, shlex

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


config = genomicsTools.scriptConfig()  #load institute-specific file server paths
netscratch = config["netscratch"]
biodata = config["biodata"]
#%%
parser = argparse.ArgumentParser()

parser.add_argument("rawReadsFolder", help="folder containing fastq.gz files for processing")
parser.add_argument("--filteredReadsFolder", "-f", help="folder for filtered reads to be put into", default="filteredreads")
parser.add_argument("--extension", "-e", help="extension of the fastq files", default=".fq.gz")
parser.add_argument("--overwrite", help="redo the filtering for files that have already been done (ie. they have a report in the folder)", action="store_true")

args = parser.parse_args()

if args.extension.startswith("."):
    extension = args.extension
else:
    extension = f".{args.extension}"
exlength = int(len(extension) * -1) - 1


filteredReadsFolder = genomicsTools.defaultFilePath(args.filteredReadsFolder, os.path.dirname(args.rawReadsFolder))

if not os.path.isdir(filteredReadsFolder):
    os.mkdir(filteredReadsFolder)

existingFiles = os.listdir(filteredReadsFolder)

with open(os.path.join(args.rawReadsFolder, "fastp_commands.txt"), "w") as commandsOut:
    
    for file in os.listdir(args.rawReadsFolder):
        if file.endswith(f"1{extension}"):
            sampleName = file.split("_")[0]
            reportFile = os.path.join(filteredReadsFolder, f"{sampleName}_fastp_report.html")
            
            if os.path.basename(reportFile) not in existingFiles or args.overwrite:
            
                inputs, outputs = [[os.path.join(folder[0], f"{sampleName}_{folder[1]}_{pair}.fq.gz") for pair in range(1, 3)] for folder in [[args.rawReadsFolder, "raw"], [filteredReadsFolder, "filt"]]]
        
                fastpCommand = genomicsTools.writeClusterCommand("fastp", (reportFile, inputs[0], inputs[1], outputs[0], outputs[1])) #returns a command to run fastp on the hpc, inserting the arugments:
                
                #bsub -q short -R "rusage[mem=1024]" -M 5012 fastp -q 20 -x --cut_front --cut_tail -W 5 -M 25 -h {0} -i {1} -I {2} -o {3} -O {4}
                
                if any([os.path.getsize(inputfile) > 1e9 for inputfile in inputs]):
                    fastpCommand = shlex.split(fastpCommand)
                    fastpCommand[2] = "normal"
                    fastpCommand[4] = '"rusage[mem=8000]"'
                    fastpCommand[6] = "10000"
                    fastpCommand = " ".join(fastpCommand)
                    
                commandsOut.write(f"{fastpCommand}\n")
                
                genomicsTools.runClusterCommand(fastpCommand, bsub = True) #runs the above command on the institute cluster
        
        
        
        
        