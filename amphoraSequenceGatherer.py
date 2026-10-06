# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 14:28:58 2022

@author: copeland
"""


import helperFunctions, os, re, logging, argparse, sys
import pandas as pd
from Bio import SeqIO, Seq

config = helperFunctions.scriptConfig() #get paths specific to our institute file server
netscratch = config["netscratch"]
biodata = config["biodata"]
logger = logging.getLogger(__name__)

logger.setLevel(logging.INFO)

seqNotFound = logging.getLogger(f"{__name__}.seqNotFound")

if __name__ == "__main__":
    interactiveArgs = {
        "strainList":[os.path.join(netscratch, "MIRO04", "MIRO04_strain_lists.xlsx"),"IS"],
        "outputPrefix":os.path.join(netscratch, "MIRO04", "MIRO04-01_IS")
         }
    parser = argparse.ArgumentParser()
    parser.add_argument("--strainList", "-s", help = "list of samples to rename the headers", nargs = "+")
    parser.add_argument("--outputPrefix", "-o", help = "folder to write results to")
    parser.add_argument("--outgroup", help = "outgroup to include in tree if not already in the strain list")
        
    args = parser.parse_args()
    
    if not hasattr(sys, "ps1"):
        logger.debug("running from command line")   
    else:
        vars(args).update(interactiveArgs)

    allMarkers = os.path.join(netscratch, "AMPHORA2", "Markers_clean") #markers already obtained from culture collection genomes using MarkerScanner.pl 
    outputFolder = args.outputPrefix

    strainList = pd.read_excel(args.strainList[0], sheet_name = args.strainList[1])["Strain_ID"].to_list()
    strainList = list(set(strainList))
    if args.outgroup:
        strainList.append(args.outgroup)
        #%%
    strainList = [entry for entry in strainList if entry != "LjR6"]

    for file in os.listdir(allMarkers):
        notFound = {item:True for item in strainList}
        allSeqs = SeqIO.parse(os.path.join(allMarkers, file), "fasta")
        with open(os.path.join(outputFolder, file), "w") as newSeqs:
            for seq in allSeqs:
                if seq.id in strainList:
                    newSeqs.write(seq.format("fasta"))
                    notFound[seq.id] = False
                    #logger.info(seq.id)
                notFoundList = [entry for entry in notFound if not entry]
                if notFoundList:
                    logger.warning(f"{', '.join(notFoundList)} not found in {file}")
    
    
    #%%
    alignCommand = helperFunctions.writeClusterCommand("MarkerAlignTrim.pl", (outputFolder,)) #writes a command line for the MarkerAlignTrim.pl inserting the arguments:
        #MarkerAlignTrim.pl -Trim -OutputFormat fasta -Directory {0}
    runAlignment = helperFunctions.runClusterCommand(alignCommand) #runs the command above on the institute cluster
    #if runAlignment.errs:
        #logger.info(runAlignment.command)
                    
    #%%                
    alignSeqs = {name:[] for name in strainList}
    
    #restoreNames = re.compile("[A-Za-z]+\d+")
    
    
    for file in os.listdir(outputFolder):
        if file.endswith("aln"):
        
            alignParse = SeqIO.parse(os.path.join(outputFolder, file), "fasta")
            align = {entry.id.split("/")[0]:entry.seq for entry in alignParse}
        
        
            alignLen = len(align[strainList[0]])
        
            logger.info(f"The alignments for {file} are {alignLen} aa long")
        
            for strainID in strainList:
                if strainID not in align:
                    logger.warning(f"{strainID} not found in {file}, using {alignLen} '-'") #if a marker is missing from a genome, replace with the correct number of "-"
                    alignSeqs[strainID].append("-"*alignLen)
                else:
                    alignSeqs[strainID].append(str(align[strainID]))
                    
               
    concatFileName = os.path.join(outputFolder, "concatenatedAlignments.fna")
    treeName = f"{os.path.splitext(concatFileName)[0]}.tre"
    #concatenated the alignments for each marker into one long sequence per strain            
    with open(concatFileName, "w") as concatFile:
        for entry in alignSeqs:
            concatFile.write(f">{entry}\n{''.join(alignSeqs[entry])}\n")
    helperFunctions.fileCreationRecord(concatFileName, vars(args), __file__)
    #%%
    makeTreeCommand = helperFunctions.writeClusterCommand("FastTreeRunner", (concatFileName, treeName)) #writes a command line for FastTree:
        #./FastTreeRunner.sh {0} {1}
    
    
    #%%
    runTree = helperFunctions.runClusterCommand(makeTreeCommand) #runs the command line above
    helperFunctions.fileCreationRecord(treeName, vars(args), __file__)        
                        
    
                    
                
            
    
