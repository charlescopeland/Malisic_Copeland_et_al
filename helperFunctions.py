# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 10:42:36 2026

@author: copeland
"""

import os
import re
import subprocess
import shlex
import json
import pandas as pd
import logging
from datetime import date

# =============================================================================
# USER CONFIGURATION
# Replace these paths with your local directories
# =============================================================================
CONFIG = {
    "netscratch": "./data/netscratch", 
    "biodata": "./data/biodata",
    "cluster": False,  # Set to True if using a Linux cluster
    "slurm": False,    # Set to True if using Slurm
    "blastbinaries": "/usr/bin", # Path to blastp, blastn, etc.
    "rdir": "/usr/bin/Rscript"    # Path to Rscript executable
}

logger = logging.getLogger(__name__)

def scriptConfig():
    """Returns the configuration dictionary."""
    return CONFIG

def defaultFilePath(fileName, defaultPath):
    """Handles path joining if the fileName is just a name or a full path."""
    if os.path.dirname(fileName):
        return fileName
    return os.path.join(defaultPath, fileName)

def parseArgList(arg):
    """Reads a list of items from a file if a file path is provided, else returns the list."""
    if isinstance(arg, str) and os.path.isfile(arg):
        with open(arg, "r") as al:
            return [entry.strip() for entry in al.readlines()]
    return arg

def readGFF(filepath, parseInfo=True, splits=[";", "="], attColumn="attribute"):
    """Parses a GFF/GTF file into a pandas DataFrame."""
    headers = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    df = pd.read_table(filepath, sep="\t", header=None, dtype=str, comment="#")
    df = df.rename(columns={i: entry for i, entry in enumerate(headers)})
    
    if parseInfo:
        def gtfParseInfo(attString):
            return {pair[0]: pair[1].replace('"', '') if len(pair) > 1 else pair[0] 
                    for pair in [pair.strip().split(splits[1]) for pair in attString.split(splits[0]) if pair]}
        df[attColumn] = df[attColumn].apply(gtfParseInfo)
    return df

def writeClusterCommand(command, options, printCommand=True, splitCommand=False):
    """
    Simulates the cluster command generation. 
    In the original, this used a template file. Here, it performs a basic 
    string format to make the scripts functional for general users.
    """
    # This is a simplified version of your template system
    # It assumes the command is a string that can be formatted with the options
    try:
        formatted_command = command.format(*options)
    except IndexError:
        formatted_command = f"{command} {' '.join(map(str, options))}"
    
    if splitCommand:
        return shlex.split(formatted_command)
    return formatted_command

def runClusterCommand(commandLine, bsub=False, returnOuts=True):
    """Executes a system command and returns the output."""
    if isinstance(commandLine, str):
        command_to_run = shlex.split(commandLine)
    else:
        command_to_run = commandLine

    # Remove cluster-specific wrappers (bsub/sbatch) for local execution
    if command_to_run[0] in ["bsub", "sbatch"]:
        # Find the actual script/command within the bsub call
        try:
            idx = command_to_run.index("./") if "./" in str(command_to_run) else 0
            command_to_run = command_to_run[idx:]
        except:
            pass

    try:
        process = subprocess.run(command_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        class Result:
            def __init__(self, stdout, stderr):
                self.outs = stdout
                self.errs = stderr
        return Result(process.stdout, process.stderr)
    except Exception as e:
        logger.error(f"Command failed: {e}")
        return type('Result', (), {'outs': '', 'errs': str(e)})

def fileCreationRecord(filePath, argsDict, script, additionalMessages=False):
    """Placeholder for the internal file history tracking system."""
    logger.debug(f"Record created for {filePath} via {script}")

def fileMetaTable(argNamespace, script=None, forExcel=True, additionalMessages=dict()):
    """Generates a metadata table for the arguments used in a run."""
    metaDict = vars(argNamespace)
    metaDict.update({"script": script, "date": date.today().strftime("%y%m%d")})
    if additionalMessages:
        metaDict.update(additionalMessages)
    return pd.DataFrame.from_dict({"arg": metaDict.keys(), "val": metaDict.values()})

def assignPlotID(plotHistoryDict, prefixLetter="M"):
    """Generates a unique ID for a plot."""
    return f"{prefixLetter}{len(plotHistoryDict) + 1:05d}"