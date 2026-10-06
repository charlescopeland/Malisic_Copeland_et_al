# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:35:02 2023

@author: copeland
"""


import helperFunctions, BioCSV, os, argparse, logging, json, sys

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


config = helperFunctions.scriptConfig()  #load institute-specific file server paths
netscratch = config["netscratch"]
biodata = config["biodata"]

interactiveArgs = {
    }

parser = argparse.ArgumentParser()

parser.add_argument("--samples", "-s", help = "samples to include in featureCounts", nargs = "+")
parser.add_argument("--reference", "-r", help = "reference")
parser.add_argument("--alignmentName","-a", help = "reference given in the alignmnent file name, if different from the reference in the .json")
parser.add_argument("--experimentFolder", "-f", help = "the folder with the experiment data")
parser.add_argument("--name", "-n", help = "an additional name for the file" )
parser.add_argument("--featureType", help = "specify a different feature type than in the reference file")

args = parser.parse_args()
if not hasattr(sys, "ps1"):
    logger.debug("running from command line")   
else:
    vars(args).update(interactiveArgs)

folder = helperFunctions.defaultFilePath(args.experimentFolder, netscratch)
featuresFolder = os.path.join(folder, "featurecounts")

reference = args.reference
with open(os.path.join(featuresFolder, "featureCounts_references.json"), "r") as rf:
    references = json.load(rf)
if reference in references:
    refInfo = references[reference]
   

#%%
if args.samples:
    suffix = ""
    if len(args.samples) == 1 and os.path.isfile((sampleFile := helperFunctions.defaultFilePath(args.samples[0], folder))):
        toCount = BioCSV.listForColumn(sampleFile, 0) #returns a list of samples from the given file
    else:
        toCount = args.samples
        
if args.name:
    suffix = args.name
    
#%%
if args.featureType:
    refInfo["featureType"] = args.featureType
    suffix = f"{suffix}_{refInfo['featureType']}"


logger.info(toCount)
logger.info(f"{len(toCount)} samples to count: {', '.join([entry for entry in toCount])}")

if args.alignmentName:
    refInBamName = args.alignmentName
else:
    refInBamName = args.reference

forFeatureCounts = [os.path.join(folder, "alignments", f"{entry}_{refInBamName}.bam") for entry in toCount]

if (noFile := [entry for entry in forFeatureCounts if not os.path.isfile(entry)]):
    logger.error(f"could not find files {noFile}")

suffix = f"_{suffix}" if suffix else ""
bsubOut = os.path.join(featuresFolder, f"bsub_out_featureCounts{suffix}.txt")
featureOut = os.path.join(featuresFolder, f"featureCounts{suffix}.tsv")
#%%
countCommand = "featureCounts" if not config["slurm"] else "runFeatureCounts"
    
featureCommand = helperFunctions.writeClusterCommand(countCommand, (refInfo["annotation"],featureOut,  refInfo["featureType"], " ".join(forFeatureCounts))) #returns a command to run featurecounts on the institute cluster:
#bsub -q normal -R "rusage[mem=5000]" -M 6000 /netscratch/dep_psl/grp_psl/CharlesSoftware/subread/bin/featureCounts -Mp --countReadPairs --primary -a {0} -A ./Genomes/chromosomeAliases.txt -o {1} -t {2} -T 4 --verbose {3}    

featureRun = helperFunctions.runClusterCommand(featureCommand) #runs the command on the institute cluster

with open(bsubOut, "w") as writeOuts:
    writeOuts.write(f"featureCounts output: {featureRun.outs}\n")
    writeOuts.write(f"featureCounts errors: {featureRun.errs}")