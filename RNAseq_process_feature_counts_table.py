# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:07:50 2023

@author: copeland
"""


import helperFunctions, os, argparse, logging, re, sys
import pandas as pd

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


config = helperFunctions.scriptConfig() #load institute-specific file server paths
netscratch = config["netscratch"]
biodata = config["biodata"]

interactiveArgs = {"inputTable":os.path.join(netscratch, "MIRO02-01", "deseq", "proc_normalized_counts.tsv"),
                   "outputTable":os.path.join(netscratch, "MIRO02-01", "kmeans", "featureCounts_sample_DESeq_normalized_sorted.tsv"),
                   #"removeSamples":os.path.join(netscratch, "VAPE15-14", "VAPE15-14_B3_samples.txt"),
                   "sampleNamePrefix" : "MM",
                   "sortOrder" : os.path.join(netscratch, "MIRO02-01", "MIRO02-01_sample_sort_order.txt")
                   }

parser = argparse.ArgumentParser()
parser.add_argument("--inputTable", "-i", help = "raw table from featureCounts")
parser.add_argument("--outputTable", "-o", help = "name for output table")
parser.add_argument("--replacementCounts", help = "count table with manual counts to replace")
parser.add_argument("--removeSamples", help = "samples to remove from the processed table", nargs = "+")
parser.add_argument("--sampleNamePrefix", help = "prefix to find sample names")
parser.add_argument("--sortOrder", help = "file with sample names in the order they should be sorted")
parser.add_argument("--fromdeseq", help = "input count table is output from the deseq script in R", action = "store_true")

args = parser.parse_args()
if not hasattr(sys, "ps1"):
    logger.debug("running from command line")   
else:
    vars(args).update(interactiveArgs)
    
folder = os.path.dirname(args.inputTable)
outputTable = helperFunctions.defaultFilePath(args.outputTable, folder)

replacementCounts = helperFunctions.defaultFilePath(args.replacementCounts, folder) if args.replacementCounts else None
sampleNamePattern = r"/(" + args.sampleNamePrefix  + r"\d{1,2})_"

with open(args.inputTable, "r") as mc:
    comments = mc.readline()
if not args.fromdeseq:
    mainDF, replaceDF = [pd.read_table(file, sep = "/t", skiprows=0, header = 1).set_index("Geneid") if file else pd.DataFrame() for file in [args.inputTable, replacementCounts] ]
else:
    mainDF, replaceDF = [pd.read_table(file, sep = " ", header = 0).set_index("Geneid") if file else pd.DataFrame() for file in [args.inputTable, replacementCounts] ]

sampleNames = re.compile(sampleNamePattern)

replaceNames = [{col : sample.group(1) for col in frame.columns if (sample := sampleNames.search(col))} for frame in [mainDF, replaceDF]]
mainDF, replaceDF = [frame.rename(columns = replaceNames[i]) for i, frame in enumerate([mainDF, replaceDF])]

for row in replaceDF.index:
    replaceData = replaceDF.loc[row]
    for x in replaceData.index:
        mainDF.loc[row, x] = replaceData.loc[x]

#%%
if args.removeSamples:
    removeSamples = [entry.strip() for entry in helperFunctions.parseArgList(args.removeSamples)]
    #treatments = pd.read_excel(os.path.join(netscratch, "MIRO02-01", "MIRO02-01_samples.xlsx"), header = 0,)
    #tDict = {entry["Sample"]:entry["Treatment"] for i, entry in treatments.iterrows()}
    #colsToRemove = [c for c in mainDF.columns if tDict.setdefault(c, True) == "frax"]
    mainDF = mainDF[[c for c in mainDF.columns if c not in removeSamples]]
    
if args.sortOrder:
    sortList = helperFunctions.parseArgList(args.sortOrder)
    sortDict = {entry:i for i, entry in enumerate(sortList)}
    allColumns = sorted(list(mainDF.columns), key = lambda x: sortDict.setdefault(x, 0))
    mainDF = mainDF[allColumns]
mainDF.to_csv(outputTable, sep = "\t", index_label = "Geneid")

helperFunctions.fileCreationRecord(outputTable, vars(args), script = __file__) #records arguments used to create the file, for my own records

        
    