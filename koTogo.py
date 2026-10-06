# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:07:50 2023

@author: copeland
"""


import helperFunctions, os, argparse, logging, math
from collections import Counter
import pandas as pd
import csv

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


config = helperFunctions.scriptConfig()
netscratch = config["netscratch"]
biodata = config["biodata"]


parser = argparse.ArgumentParser()
parser.add_argument("--organism", "-o", help = "organisms to get data for")
parser.add_argument("--folder", "-f", help = "folder to save the results")
parser.add_argument("--requeryKO", help = "get the mapping of genes and KO terms again from the KEGG website", action="store_true")
parser.add_argument("--duplicateKOs", help = '''how to deal with KO terms that apply to multiple genes with different GO terms. 
                    options are 
                    'all': KO term is assigned all GO terms associated with any gene it covers,
                    'vote': KO term is assigned GO terms associated with at least 50% of the genes it covers,
                    'shared': KO term is assigned only the common GO terms shared with every genes it covers''', default = "all")

args = parser.parse_args(["-o", "pae", "-f", os.path.join(netscratch, "MIRO02-01", "functions"), "--duplicateKOs", "low"])


import requests
import time

ORG = args.organism#"pae"  # KEGG organism code for Pseudomonas aeruginosa (PAO1)
folder = args.folder
def padLocus(locus):
    return f"{'_'.join(locus.split('_')[:-1])}_{locus.split('_')[-1].zfill(4)}"

def get(url, pause=1):
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    time.sleep(pause)  # be nice to KEGG servers
    return r.text

def parse_link_tsv(text):
    # Parses KEGG "link" TSV: left<TAB>right per line
    pairs = []
    for line in text.strip().splitlines():
        if not line.strip():
            continue
        a, b = line.split("\t")
        pairs.append((a, b))
    return pairs

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]
        
def majority_items(series, frac=0.5):
    # Count each item at most once per row
    lists = [set(x) for x in series.dropna()]
    if not lists:
        return []
    n = len(lists)
    thr = math.ceil(frac * n)   # "at least 50%" -> ceil
    cnt = Counter()
    for s in lists:
        cnt.update(s)
    return [item for item, c in cnt.items() if c >= thr]



# 1) Get gene→KO for the organism
genetokofile = os.path.join(folder, f"{ORG}_gene_to_ko.tsv")
if args.requeryKO:
    print(f"Downloading gene→KO for organism {ORG})...")
    gene_ko_txt = get(f"https://rest.kegg.jp/link/ko/{ORG}")
    gene_ko_pairs = parse_link_tsv(gene_ko_txt)  # e.g., ("pae:PA0001", "ko:Kxxxxx")
    
    # Normalize
    gene_ko = []
    for g, k in gene_ko_pairs:
        gene = g.split(":", 1)[1]  # PAxxxx
        ko = k.split(":", 1)[1]    # Kxxxxx
        gene_ko.append((gene, ko))
    
    # Save KO to gene name (PAxxxx)
    
    with open(genetokofile, "w") as f:
        f.write("gene\tko\n")
        for gene, ko in sorted(set(gene_ko)):
            f.write(f"{gene}\t{ko}\n")
    print(f"gene_to_ko.tsv: {len(set(gene_ko))} pairs")

# PA to GO

geneToGo = pd.read_csv(os.path.join(folder, "pseudocap.gaf"), sep="\t", comment="!", header = None)

geneToGo = geneToGo.dropna(subset = 10).copy()
geneToGo["PA"] = geneToGo[10].apply(lambda x : x.split("|")[-1])

godict = {}

for i, row in geneToGo[["PA", 4]].iterrows(): #BiologicalProcesses.iterrows():
    godict.setdefault(row["PA"], []).append(row[4])

godict = {key : list(set(val)) for key, val in godict.items()}


KOtoPA = pd.read_table(genetokofile, sep = "\t")

KOtoPA["GO"] = KOtoPA["gene"].apply(lambda x: godict.setdefault(x, []))
thresholds = {"vote":0.5, "shared":1, "low":0.02}
if (duplicateKOs := args.duplicateKOs.casefold()).startswith("a"):
    KOdict = {}
    for i, row in KOtoPA.iterrows():
        KOdict.setdefault(row["ko"], [])
        for entry in row["GO"]:
            if entry not in KOdict[row["ko"]]:
                KOdict[row["ko"]].append(entry)
elif duplicateKOs in thresholds:
    KOdict = KOtoPA.groupby('ko')['GO'].apply(majority_items, frac=thresholds[duplicateKOs]).to_dict()
else:
    raise NameError("duplicateKOs argument {args.duplicateKOs} not recognized")            
#KOdict = dict(zip(KOtoPA["ko"], KOtoPA["GO"]))
#%%
# IT466 genes to GO, via KO > PA

IT466 = helperFunctions.readGFF(os.path.join(netscratch, "Genomes", "IT466.gff"))

IT466["KO"] = IT466["attribute"].apply(lambda x: x.get("KO"))

IT466["GO"] = IT466["KO"].apply(lambda x: KOdict.get(x))

IT466["GOstr"] = IT466["GO"].apply(lambda x: ", ".join(x) if x else "")


IT466["gene"] = IT466["attribute"].apply(dict.get, args = ("ID", ))
IT466["gene"] = IT466["gene"].apply(padLocus)


IT466[["gene", "GOstr"]].to_csv(os.path.join(folder, f"IT466_gene_to_go_from_aeruginosa_{args.duplicateKOs}.map"), sep = "\t", quoting = csv.QUOTE_NONE, header = False, index = False)

# parse and add GO terms from interproscan analysis
gff = helperFunctions.readGFF(os.path.join(folder, "IT466_interpro_pfam.gff3"))

gff["gene"] = gff["chromosome"].apply(padLocus)
def parseOntology(attdict):
    if (GOstring := attdict.get("Ontology_term")):
        return [term.replace('"', '').strip() for term in GOstring.split(",")]
    
gff["GO"] = gff["attribute"].apply(parseOntology)

interproGO = {}
for i, row in gff.iterrows():
    if row["GO"]:
        [interproGO.setdefault(row["gene"], []).append(entry) for entry in row["GO"]]
    
dictForFrame = {}
dictForFrame["gene"] = list(interproGO.keys())
dictForFrame["GOstr"] = [', '.join(entry) if entry else "" for entry in interproGO.values()]

mappingFrame = pd.DataFrame(data = dictForFrame)
    
if False:
    mappingFrame[["gene", "GOstr"]].to_csv(os.path.join(folder, "IT466_gene_to_go_interpro.map"), sep = "\t", quoting = csv.QUOTE_NONE, header = False, index = False)

combinedFrame = IT466[["gene", "GO"]].copy()
combinedFrame["GO"] = combinedFrame["GO"].apply(lambda x: x if isinstance(x, list) else [])
combinedFrame["combined"] = combinedFrame.apply(lambda x: list(set(x["GO"] + interproGO.setdefault(x["gene"], []))), axis = 1)
combinedFrame["GOstr"] = combinedFrame["combined"].apply(lambda x: ", ".join(x) if x else "")

combinedFrame[["gene", "GOstr"]].to_csv(os.path.join(folder, f"IT466_gene_to_go_{args.duplicateKOs}.map"), sep = "\t", quoting = csv.QUOTE_NONE, header = False, index = False)

