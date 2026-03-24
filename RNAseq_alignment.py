# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:03:30 2023

@author: copeland
"""


import genomicsTools, os, logging, argparse, BioCSV, re, time
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


config = genomicsTools.scriptConfig()  #load institute-specific file server paths
netscratch = config["netscratch"]
biodata = config["biodata"]

parser = argparse.ArgumentParser()

parser.add_argument("--reference", "-r", help="prefix of hisat2 index to align reads against")
parser.add_argument("--sample", "-s", help="sample to align")
parser.add_argument("--file", "-f", help="csv file containing a samples and corresponding reference to align against. If given, ignores other sample and reference arguments")
parser.add_argument("--outputSuffix", "-o" , help="extra suffix to add to output files")
parser.add_argument("--directory", "-d", help = "directory with reads to align")
parser.add_argument("--noSplicedAlignment", help = "use the --no-spliced-alignment option in hisat2", action = "store_true")

args = parser.parse_args()

minutes = 60

if args.file:
    toAlign = BioCSV.listTotal(args.file) #returns list of files from file
else:
    toAlign = [[args.sample,  args.reference]]
    
if not os.path.isdir(alignmentFolder := os.path.join(args.directory, "alignments")):
    os.mkdir(alignmentFolder)

if args.outputSuffix:
    suff = f"_{args.outputSuffix}"
else:
    suff = ""

if len(toAlign) <= 4:
    sleepTime = 0.1*minutes
else:
    sleepTime = 1 * minutes
    
for sample, reference in toAlign:
    fastqs = [os.path.join(args.directory, "filteredreads", f"{sample}_filt_{mate}.fq.gz") for mate in [1,2]]
    
    refName = re.match("(.*)_hisat_ref", os.path.basename(reference))
    if refName:
        refName = refName.group(1)
    else:
        refName =  os.path.basename(toAlign[sample])
    outputPref = os.path.join(alignmentFolder, f"{sample}_{refName}{suff}") 
    commandString = "hisat2_alignment.sh" if not args.noSplicedAlignment else "hisat2_alignment_nosplice.sh"
    alignmentCommand = genomicsTools.writeClusterCommand(commandString, (outputPref, reference) + tuple(fastqs)) #returns command to run on institute cluster, inserting given arguments:
        #bsub -q normal -R "rusage[mem=4000]" -M 8000 ./hisat2_alignment.sh -p {0} {1} {2} {3}
        #or
        #bsub -q normal -R "rusage[mem=4000]" -M 8000 ./hisat2_alignment_nosplice.sh -p {0} {1} {2} {3}
        
    logger.debug(f"alignment for {sample} running with: {alignmentCommand}")
    runAlignment = genomicsTools.runClusterCommand(alignmentCommand, bsub = f"bsub_out_{sample}_{refName}{suff}_alignment.txt", returnOuts = False) #runs the above command on the institute cluster
    
    time.sleep(sleepTime)
    
    
    logger.debug(runAlignment.outs)
    logger.debug(runAlignment.errs)