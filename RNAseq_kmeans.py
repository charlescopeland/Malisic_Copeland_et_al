# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:07:50 2023

@author: copeland
"""


import genomicsTools, os, argparse, logging, sys, json
import scipy.cluster.vq as spclus
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import heatmap
import figureImproveR, pythonPlots
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
config = genomicsTools.scriptConfig()
netscratch = config["netscratch"]
biodata = config["biodata"]
cm = 1 / 2.5
def scale_row(row):
    mean = row.mean()
    std = row.std()
    return (row - mean) / std
plotHistoryFile = os.path.join(biodata, "R_scripts", "pythonScripts", "heatmap_history.json")
with open(plotHistoryFile, "r") as jf:
    previousPlots = json.load(jf)
    
if __name__ == "__main__":
    interactiveArgs = {"countTableFile":os.path.join(netscratch, "MIRO02-01", "kmeans", "featureCounts_sample_DESeq_normalized_sorted_coumarin_top.tsv"),
                       "geneList":os.path.join(netscratch, "MIRO02-01", "deseq", "wt_specific_gene_list_unpadded.txt"),
                       "outputPrefix":"MIRO02-01_wt-specific_genes_hierarch",
                       "treatments" : [os.path.join(netscratch, "MIRO02-01", "MIRO02-01_samples.xlsx" ), "samples", "Sample"],
                       "chosenClusters":1,
                       "clusteringType":"hierarchical",
                       "fromdeseq":True,
                       "purpose":"publication",
                       "figureSize": [8, 8],
                       "gapSpacing" : 0,
                       "outputClusters":False,
                       "aspect":1
         }
    
    if interactiveArgs.get("update"):
        interactiveArgs = {"update":interactiveArgs["update"],
                                                      }
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--countTableFile", "-f", help = "table of counts to cluster")
    parser.add_argument("--geneList", "-g", help = "file with list of genes to include")
    parser.add_argument("--outputPrefix", "-o", help = "prefix for output files")
    parser.add_argument("--treatments", help = "table with treatments for the samples, including file path, sheet name, and optional name of column with the samples", nargs = "+")
    parser.add_argument("--chosenClusters", help = "number of clusters to use", type = int)
    parser.add_argument("--clusteringType", help = "use k-means or hierarchical clustering", default = "kmeans")
    parser.add_argument("--fromdeseq", help = "count file is normalized by DESeq2, and lacks the first columns from featureCounts", action = "store_true")
    parser.add_argument("--purpose", help = "presentation or publication")
    parser.add_argument("--figureSize", help = "size of the heatmap", nargs = "+")
    parser.add_argument("--aspect", help = "aspect ratio for the heatmap")
    parser.add_argument("--replicate", help = "replicate an existing heatmap")
    parser.add_argument("--update", help = "create a new figure based on an existing one, changing only a few options")
    parser.add_argument("--outputClusters", help = "output individual text files listing the genes in each cluster", action = "store_true"),
    parser.add_argument("--gapSpacing", help = "number of white rows to include as gaps between clusters", default = 25)
    args = parser.parse_args()
    
    if not hasattr(sys, "ps1"):
        logger.debug("running from command line")   
    else:
        vars(args).update(interactiveArgs)
    vars(args)["outputPrefix"] = f"{vars(args)['outputPrefix']}_k{args.chosenClusters}"  
    newPlotid = genomicsTools.assignPlotID(previousPlots, "M")
    if (replicate := args.replicate):
        try: 
            logger.info(f"Replicating figure {args.replicate}")
            args = argparse.Namespace(**previousPlots[args.replicate])
        except KeyError as noRep:
            raise KeyError(f"figure id {noRep} not found")
        vars(args).update({"plotID":args.replicate})
    elif args.update:
        try:
            logger.info(f"Basing figure on {args.update}")
            previousArgs = previousPlots[args.update]
            if not "outputPrefix" in interactiveArgs:
                newSuffix = f"{previousArgs['outputPrefix']}_{newPlotid}"
                logger.warning(f"no new figure suffix given. saving figure with suffix {newSuffix}")
                interactiveArgs["figureSuffix"] = newSuffix
            previousArgs.update(interactiveArgs)
            
            args = argparse.Namespace(**previousArgs)
            vars(args).update({"plotID":newPlotid}) 
        except KeyError:
            raise KeyError(f"figure id {args.replicate} not found")    
    
 #%%       
    
    counts = pd.read_csv(args.countTableFile, sep="\t", header=0, comment="#", index_col=0)

    plt.style.use(os.path.join(biodata, "R_scripts", "pythonScripts", f"mplstyle.{args.purpose}.txt"))
    if not args.fromdeseq:
        counts = counts.iloc[:, 5:]
    degs = genomicsTools.parseArgList(args.geneList)
    args.outputPrefix = genomicsTools.defaultFilePath(args.outputPrefix, os.path.dirname(args.countTableFile))
    ## the following is modified from chatGPT
    # 1) Subset to significant DEGs
    if degs:
        use_counts = counts.loc[counts.index.intersection(degs)].copy()
    else:
        use_counts = counts.copy()
    
    cwm = matplotlib.colors.LinearSegmentedColormap.from_list(
        'cyan_white_magenta',
        [(0.0, '#00BFFF'),   # cyan for low (negative)
         (0.5, '#FFFFFF'),   # white at center (zero)
         (1.0, '#8800aa')],  # magenta for high (positive)
        N=256
    ).with_extremes(bad='white')  # NaNs show as white (for gap rows, etc.)
    v = 2
    norm = matplotlib.colors.TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)
    # Optional: drop genes with extremely low total counts (defensive)
    use_counts = use_counts[use_counts.sum(axis=1) > 0]
    if False:
        # 2) Normalize to CPM and log-transform
        def log2_cpm(df):
            lib_sizes = df.sum(axis=1)  # careful: rows are genes; we want library sizes per sample (columns)
            # Correct: sum across rows for each column (per-sample library size)
            lib_sizes = df.sum(axis=0)
            cpm = df.div(lib_sizes, axis=1) * 1e6
            return np.log2(cpm + 1)
        
        logcpm = log2_cpm(use_counts)
        
        # 3) Per-gene z-score to focus on patterns across samples
        gene_means = logcpm.mean(axis=1)
        gene_stds = logcpm.std(axis=1)
        Z = logcpm.sub(gene_means, axis=0).div(gene_stds.replace(0, np.nan), axis=0).dropna()
        
        X = Z.values.astype(float)  # genes as rows (observations), samples as columns (features)
    else:
        Z = np.log2(use_counts + 1).apply(scale_row, axis = 1)
        X = Z.values.astype(float)
    if args.clusteringType.casefold().startswith("k"):
        # 4) Pick k (elbow on distortion)
        ks = range(2, 11)
        distortions = []
        for k in ks:
            # kmeans returns (centroids, distortion); set more iterations for stability
            _, dist = spclus.kmeans(X, k, iter=50, thresh=1e-5)
            distortions.append(dist)
        
        plt.figure(figsize=(4,3))
        plt.plot(ks, distortions, 'o-')
        plt.xlabel('k')
        plt.ylabel('Distortion (within-cluster SSE proxy)')
        plt.title('Elbow plot')
        plt.tight_layout()
        plt.show()
        #%%
        # Choose k based on the elbow (replace with your chosen k)
        k = args.chosenClusters
    
        # 5) Final k-means clustering with kmeans2 (gets labels)
        centroids, labels = spclus.kmeans2(X, k, minit='++', iter=100, seed=0)
        
        # 6) Results as DataFrames
        clusters = pd.Series(labels, index=Z.index, name='cluster')
        centroids_df = pd.DataFrame(centroids, columns=Z.columns)  # cluster profiles in z-score space
        
            
        # 7) Example: attach cluster labels to your DEGs
        deg_clustered = Z.copy()
        deg_clustered['cluster'] = clusters
        
        # 8) (Optional) visualize cluster centroids
        plt.figure(figsize=(6,3))
        for i in range(k):
            plt.plot(centroids_df.columns, centroids_df.iloc[i], alpha=0.6, label=f'Cluster {i}')
        plt.legend(ncol=2, fontsize=8)
        plt.xlabel('Samples')
        plt.ylabel('Z-score')
        plt.title('Cluster centroids (z-scored)')
        plt.tight_layout()
        
        
        with pd.ExcelWriter(f"{args.outputPrefix}.xlsx", mode = "w") as ew:
            centroids_df.to_excel(ew, sheet_name="centroids")
            deg_clustered.to_excel(ew, sheet_name="clusters")
            #%%
            if args.treatments:
                treatments = pd.read_excel(args.treatments[0], sheet_name=args.treatments[1])
                if len(args.treatments) > 2:
                    treatments = treatments.set_index(args.treatments[2])
                else:
                    treatments.set_index(0)
                treatMerge = pd.merge(centroids_df.transpose().rename(columns={entry:f"cluster_{entry}" for entry in range(0, args.chosenClusters)}), treatments, left_index = True, right_index = True)
                treatMerge.to_excel(ew, sheet_name="toPlot")
            genomicsTools.fileMetaTable(args, script = __file__).to_excel(ew, sheet_name = "meta")
            genomicsTools.fileCreationRecord(f"{args.outputPrefix}.xlsx", vars(args), script = __file__)
    elif args.clusteringType.casefold().startswith("h"):
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    # 1) Hierarchical clustering of genes (rows)
    # Use correlation distance (1 - Pearson r) and average linkage (common for expression data)
        dist_genes = pdist(Z.values, metric='correlation')
        link_genes = linkage(dist_genes, method='ward', optimal_ordering=True)
        #%%
        if args.figureSize:
            figWidth, figHeight = [float(entry) * cm for entry in args.figureSize]
        else:
            figWidth = 8
            figHeight = 8
        dendFigure = plt.figure(figsize=(figWidth, figHeight))
        dendro_genes = dendrogram(
            link_genes,
            labels=Z.index.to_list(),
            leaf_rotation=0,
            leaf_font_size=6,
            orientation="left"
        )
        dendFigure.suptitle('Gene dendrogram (correlation distance, average linkage)')
        dendFigure.supylabel('1 - Pearson correlation')
        dendFigure.tight_layout()
        #dendFigure.show()
 #%%       
        # 2) (Optional) Cut the dendrogram into k flat clusters
        k = args.chosenClusters
        gene_labels = fcluster(link_genes, t=k, criterion='maxclust')
        gene_clusters = pd.Series(gene_labels, index=Z.index, name='cluster')
        # dendFigure = plt.figure(figsize=(6, 4))
        
        # dendro_samples = dendrogram(
        #     link_samples,
        #     labels=Z.columns.to_list(),
        #     leaf_rotation=90
        # )
        # plt.title('Sample dendrogram (correlation distance, average linkage)')
        # plt.tight_layout()
        # plt.show()
        
        # 4) (Optional) Heatmap ordered by the dendrogram leaves (no seaborn)
        # Get leaf orders without replotting for samples
        #sample_leaves = dendrogram(link_samples, no_plot=True)['leaves']
        gene_order = dendro_genes['ivl']  # labels already ordered by the gene dendrogram
        gene_leaves = dendro_genes['leaves']        # indices in the order of the dendrogram
        Z_ord = Z.loc[gene_order, :]
        clusters = pd.merge(Z_ord, gene_clusters, left_index=True, right_index=True)
        #%%
        plotGap = args.gapSpacing
        blocks = []
        for i, c in enumerate(clusters["cluster"].unique()):
            g = clusters.loc[clusters["cluster"] == c].copy()
            blocks.append(g)
            if i < len(clusters) - 1:
                gap_index = [f'gap_{i}_{j}' for j in range(plotGap)]
                blocks.append(pd.DataFrame(np.nan, index=gap_index, columns=Z.columns))               
        Z_gap = pd.concat(blocks, axis=0)
        clusterAnnotation = Z_gap.loc[:, ["cluster"]]
        Z_gap = Z_gap[[entry for entry in Z_gap.columns if entry != "cluster"]]
        
        
        #%%
        fig = plt.figure(figsize=(8, 8), constrained_layout=True)
        gs = fig.add_gridspec(1, 2, width_ratios=[20, 1], wspace=0.05)
        ax_hm = fig.add_subplot(gs[0, 0])
        im = heatmap.heatmap(Z_gap, row_labels=False, col_labels=Z_gap.columns, ax = ax_hm, cbarlabel="Scaled Counts", aspect = args.aspect, grid=False, cmap = cwm, norm=norm, interpolation="none",
                             cbar_kw = {"shrink":0.1, "aspect" : 10})
        ax_an = fig.add_subplot(gs[0, 1], sharey=ax_hm)
        cluster_cmap = plt.get_cmap('tab20', k)  # or define your own palette
        cluster_cmap = matplotlib.colors.ListedColormap([cluster_cmap(i) for i in range(k)]).with_extremes(bad='white')
        annim = heatmap.heatmap(clusterAnnotation, row_labels = False, col_labels = ["cluster"], cmap = cluster_cmap, aspect = "auto", grid = False)
        [axis.get_yaxis().set_visible(False) for axis in [ax_hm, ax_an]]
        
        
        
        #%%
        #im = plt.imshow(Z_gap, aspect='auto', cmap=cwm, norm=norm, interpolation="none")# vmin=-2.0, vmax=2.0)
        #im = plt.imshow(Z_ord, aspect='auto', cmap='seismic', vmin=-2.0, vmax=2.0, interpolation="none")
        #plt.colorbar(im, label='scaled counts', shrink = 0.2, aspect = 10)
        #plt.xlabel('Samples')
        #fig.subplots().set_ylabel('Genes (hierarchically ordered)')
        #fig.suptitle('DEG expression heatmap (z-scored)')
        #fig.tight_layout()
        #plt.show()
        #%%
        hmFileName, dendFileName = [f"{args.outputPrefix}_{figType}.svg" for figType in ["heatmap", "dendrogram"]]
        fig.savefig(hmFileName)
        dendFigure.savefig(dendFileName)
        replacement = figureImproveR.commonTextReplacement()
        [figureImproveR.improveFigure(file, resize = False, textHeight=True, textReplace = replacement)  for file in [hmFileName, dendFileName]]
        
        metadata = [f"Counts from {args.countTableFile}",
                    f"Genes from {args.geneList}",
                    f"Treatments from {args.treatments}",
                    f"{args.clusteringType} clustering with k = {args.chosenClusters}",
                    f"figure size: {figWidth} x {figHeight}",
                    f"heatmap file : {hmFileName}",
                    f"dendrogram file: {dendFileName}",
                    f"plot ID: {newPlotid}"]         
        [pythonPlots.addMetadata(file, metadata, startingPosition = [0, figHeight + 2*cm]) for file in [hmFileName, dendFileName]]
        with pd.ExcelWriter(f"{args.outputPrefix}.xlsx", mode = "w") as ew:
            Z_ord.to_excel(ew, sheet_name="gene_order")
            
            clusters.to_excel(ew, sheet_name="clusters")
            centroids = clusters.groupby("cluster").mean()
            centroids.to_excel(ew, sheet_name = "centroids")
            #%%
            if args.treatments:
                treatments = pd.read_excel(args.treatments[0], sheet_name=args.treatments[1])
                if len(args.treatments) > 2:
                    treatments = treatments.set_index(args.treatments[2])
                else:
                    treatments.set_index(0)
                treatMerge = pd.merge(centroids.transpose().rename(columns={entry:f"cluster_{entry}" for entry in range(0, args.chosenClusters + 2)}), treatments, left_index = True, right_index = True)
                treatMerge.to_excel(ew, sheet_name="toPlot")
            genomicsTools.fileMetaTable(args, script = __file__).to_excel(ew, sheet_name = "meta")
            genomicsTools.fileCreationRecord(f"{args.outputPrefix}.xlsx", vars(args), script = __file__)
            
        if not replicate:
            previousPlots.update({newPlotid: vars(args)})
            with open(plotHistoryFile, "w") as phf:
                json.dump(previousPlots, phf, indent = 4)
        

 #%%               
    if args.outputClusters:
        geneFiles = []
        for g in clusters["cluster"].unique():
            genes = clusters.loc[clusters["cluster"] == g].index.to_list()
            geneFileName = f"{args.outputPrefix}_cluster_{g}_genes.txt"
            with open(geneFileName, "w") as gf:
                gf.write("\n".join([str(gene) for gene in genes]))
            geneFiles.append(geneFileName)
        fofName = f"{args.outputPrefix}_cluster_genes.fof"
        with open(fofName, "w") as fof:
            fof.write("\n".join(geneFiles))
        fofArgs = vars(args).copy()
        fofArgs.update({"files included": geneFiles})
        genomicsTools.fileCreationRecord(fofName, fofArgs, script = __file__ )
            
       