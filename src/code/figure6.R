library(scater)
library(scRNAseq)
library(Rtsne)
require(fastcluster)
library(SingleCellExperiment)
set.seed(5252)

# Data
example_sce <- ZeiselBrainData()
example_sce
example_sce <- addPerCellQC(example_sce, 
                            subsets=list(Mito=grep("mt-", rownames(example_sce))))
plotColData(example_sce, x = "sum", y="detected", colour_by="tissue") 

example_sce <- logNormCounts(example_sce)

plotExpression(example_sce, rownames(example_sce)[1:6], x = "level1class")

# ACP
example_sce <- runPCA(example_sce)
plotPCA(example_sce, colour_by="Snap25")
example_sce$assays
colData(example_sce)


#test
#pvals[b] <-  test_hier_clusters_exact(X, link="average", K=3, k1=1, k2=2, hcl=hcl)$pval
